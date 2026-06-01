"""Multi-clone aggregation of samtools mpileup outputs to compute a Levanon-style DEI.

After Path A: N=4 BE4 clones (BE4_clone1..4) + M=3 nCas9 controls + Parent_WGS.
Multi-clone aggregation suppresses clone-specific somatic noise that contaminated the single-clone run.

Per gene:
  DEI_tpc = sum_{treated_clones} edited_reads / sum_{treated_clones} total_reads
          - sum_{control_clones}   edited_reads / sum_{control_clones}   total_reads
  (Levanon-style pooled rate, treated minus matched control)

Shared whitelist: drop position iff ANY control sample shows VAF > 0.05.
"""
import sys, os, pandas as pd, numpy as np
from collections import defaultdict

OUTDIR = '/mnt/data/dna_features/index'
GOLD = '/mnt/data/gold_index'
TPC_BED = f'{OUTDIR}/cds_tpc.bed'
NPC_BED = f'{OUTDIR}/cds_npc.bed'
PARENT_VAF_CAP = 0.05
MIN_POS_PER_GENE = 20

# Treated and control sample sets
TREATED = ['BE4_clone1', 'BE4_clone2', 'BE4_clone3', 'BE4_clone4']
CONTROL = ['nCas9_clone1', 'nCas9_clone2', 'nCas9_clone3', 'Parent_WGS']

# Load annotation
ann = {}
for bedp, mot in [(TPC_BED, 'tpc'), (NPC_BED, 'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom, s, e, gene, sc, strand = f
            ann[(chrom, int(s) + 1)] = (gene, mot, strand)
print(f"loaded {len(ann):,} annotated positions", file=sys.stderr)

def parse_bases(bases, ref):
    out = {'A':0,'C':0,'G':0,'T':0}
    R = ref.upper(); i = 0; n = len(bases)
    while i < n:
        c = bases[i]
        if c == '^': i += 2; continue
        if c == '$': i += 1; continue
        if c in '+-':
            j = i + 1
            while j < n and bases[j].isdigit(): j += 1
            try: i = j + int(bases[i+1:j])
            except ValueError: i += 1
            continue
        if c == '*': i += 1; continue
        if c in '.,':
            if R in out: out[R] += 1
        elif c.upper() in 'ACGT':
            out[c.upper()] += 1
        i += 1
    return out

def parse_pileup(path, label):
    """Return dict (chrom,pos) -> {tot, main, noise, v_main, strand, motif, gene}"""
    out = {}; n_rows = 0
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5: continue
            n_rows += 1
            chrom = f[0]; pos = int(f[1]); ref = f[2].upper()
            depth = int(f[3])
            if depth < 1: continue
            key = (chrom, pos)
            if key not in ann: continue
            gene, motif, strand = ann[key]
            ac = parse_bases(f[4], ref)
            tot = sum(ac.values())
            if tot < 1: continue
            if strand == '+':
                if ref != 'C': continue
                main = ac['T']; noise = (ac['A'] + ac['G']) / 2.0
            else:
                if ref != 'G': continue
                main = ac['A']; noise = (ac['C'] + ac['T']) / 2.0
            out[key] = (tot, main, noise, main/tot, motif, gene)
            if n_rows % 2_000_000 == 0:
                print(f"  {label}: {n_rows:,} rows -> {len(out):,} kept", file=sys.stderr)
    print(f"  {label}: total {n_rows:,} kept {len(out):,}", file=sys.stderr)
    return out

# Load all samples (only those whose pileup file exists)
samples = {}
for s in TREATED + CONTROL:
    pile = f'{GOLD}/{s}_samtools.mpileup'
    if not os.path.exists(pile):
        print(f"  missing: {pile}", file=sys.stderr); continue
    print(f"parsing {s}...", file=sys.stderr)
    samples[s] = parse_pileup(pile, s)

treated_samples = [s for s in TREATED if s in samples]
control_samples = [s for s in CONTROL if s in samples]
print(f"\ntreated samples loaded: {treated_samples}", file=sys.stderr)
print(f"control samples loaded: {control_samples}", file=sys.stderr)

# All positions present in ALL samples (intersection)
common = set.intersection(*[set(samples[s].keys()) for s in treated_samples + control_samples])
print(f"positions in ALL samples: {len(common):,}", file=sys.stderr)

# Whitelist: drop position if ANY control sample shows VAF > cap
whitelist = set(); dropped = 0
for key in common:
    if any(samples[c][key][3] > PARENT_VAF_CAP for c in control_samples):
        dropped += 1; continue
    whitelist.add(key)
print(f"after control-VAF filter: dropped={dropped:,} ({100*dropped/max(len(common),1):.2f}%); kept={len(whitelist):,}", file=sys.stderr)

# Per-gene multi-clone aggregation
# treated_main = sum over (positions in whitelist) of sum over treated_clones of main
# treated_tot  = sum over (positions in whitelist) of sum over treated_clones of tot
agg = defaultdict(lambda: {
    'tT_pos':0, 'tT_main':0, 'tT_tot':0, 'tT_noise':0,
    'tC_pos':0, 'tC_main':0, 'tC_tot':0, 'tC_noise':0,
    'nT_pos':0, 'nT_main':0, 'nT_tot':0, 'nT_noise':0,
    'nC_pos':0, 'nC_main':0, 'nC_tot':0, 'nC_noise':0,
})
for key in whitelist:
    motif, gene = samples[treated_samples[0]][key][4], samples[treated_samples[0]][key][5]
    p = 't' if motif == 'tpc' else 'n'
    d = agg[gene]
    d[f'{p}T_pos'] += 1  # increments once per position (not per clone)
    d[f'{p}C_pos'] += 1
    for s in treated_samples:
        tot, main, noise, _, _, _ = samples[s][key]
        d[f'{p}T_tot'] += tot; d[f'{p}T_main'] += main; d[f'{p}T_noise'] += noise
    for s in control_samples:
        tot, main, noise, _, _, _ = samples[s][key]
        d[f'{p}C_tot'] += tot; d[f'{p}C_main'] += main; d[f'{p}C_noise'] += noise

rows = []
for gene, d in agg.items():
    row = {'gene': gene}
    for label in ('tT','tC','nT','nC'):
        denom = max(d[f'{label}_tot'], 1)
        row[f'{label}_n_pos'] = d[f'{label}_pos']
        row[f'{label}_depth'] = d[f'{label}_tot']
        row[f'{label}_EI_main'] = d[f'{label}_main'] / denom
        row[f'{label}_EI_noise'] = d[f'{label}_noise'] / denom
    rows.append(row)
df = pd.DataFrame(rows)
df['DEI_tpc'] = df.tT_EI_main - df.tC_EI_main
df['DEI_npc'] = df.nT_EI_main - df.nC_EI_main
df['noise_gap_tpc'] = df.tT_EI_main - df.tT_EI_noise
df['noise_gap_npc'] = df.nT_EI_main - df.nT_EI_noise

df_full = df.copy()
df = df[(df.tT_n_pos >= MIN_POS_PER_GENE) & (df.tC_n_pos >= MIN_POS_PER_GENE)].copy()
df.to_parquet(f'{OUTDIR}/per_gene_dei_multiclone.parquet')

def s(a): return f'n={len(a)} mean={a.mean():.5g} med={a.median():.5g} q90={a.quantile(0.9):.5g}'
qc = [
    f"=== Multi-clone DEI QC ===",
    f"treated samples: {treated_samples}",
    f"control samples: {control_samples}",
    f"PARENT_VAF_CAP={PARENT_VAF_CAP} MIN_POS_PER_GENE={MIN_POS_PER_GENE}",
    f"positions in all samples: {len(common):,} | dropped (any control VAF > {PARENT_VAF_CAP}): {dropped:,} | whitelist: {len(whitelist):,}",
    f"genes after coverage filter: {len(df):,} of {len(df_full):,}",
    f"",
    f"--- Per-gene rates (multi-clone pooled) ---",
    f"Treated EI_main TpC: {s(df.tT_EI_main)}",
    f"Control EI_main TpC: {s(df.tC_EI_main)}",
    f"Treated EI_main nonTpC: {s(df.nT_EI_main)}",
    f"Control EI_main nonTpC: {s(df.nC_EI_main)}",
    f"",
    f"--- DEI (HEADLINE) ---",
    f"DEI_tpc: {s(df.DEI_tpc)}",
    f"DEI_npc (motif-NEG): {s(df.DEI_npc)}",
    f"",
    f"--- Pre-registered gates ---",
    f"[a] DEI_tpc mean > 0: {df.DEI_tpc.mean() > 0}",
    f"[b] Mean(DEI_tpc)/Mean(DEI_npc): {df.DEI_tpc.mean() / max(df.DEI_npc.mean(), 1e-12):.2f} (need > 2x convincing / > 1.3x weak)",
    f"[c] Frac genes DEI_tpc > DEI_npc: {(df.DEI_tpc > df.DEI_npc).mean():.3f} (need > 0.70)",
    f"[d] Frac genes Treated TpC > Control TpC: {(df.tT_EI_main > df.tC_EI_main).mean():.3f}",
    f"",
    f"--- Top-100 motif specificity (the strongest single test) ---",
    f"Top-100 DEI_tpc mean: {df.nlargest(100,'DEI_tpc').DEI_tpc.mean():.5g}",
    f"Top-100 DEI_npc mean: {df.nlargest(100,'DEI_tpc').DEI_npc.mean():.5g}",
    f"Top-100 ratio: {df.nlargest(100,'DEI_tpc').DEI_tpc.mean() / max(df.nlargest(100,'DEI_tpc').DEI_npc.mean(), 1e-12):.2f}x",
]
text = '\n'.join(qc)
print(text)
with open(f'{OUTDIR}/qc_summary_multiclone.txt', 'w') as fh: fh.write(text + '\n')
print(f"\nwrote {OUTDIR}/per_gene_dei_multiclone.parquet", file=sys.stderr)
