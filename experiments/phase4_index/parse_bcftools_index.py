"""Parse pre-existing /mnt/data/gold_index/{BE4,Parent}_chr*.tsv (bcftools mpileup -a AD output)
into per-position counts at TpC + non-TpC sites, then compute per-gene DNA Editing Index.

QA-driven design (post-review):
 - QA-H3 fix: SHARED VAF whitelist across treated AND control. A position with VAF > VAF_CAP
   in EITHER sample is dropped from BOTH samples. Prevents asymmetric germline-SNV bias.
 - QA-H2 fix: NO double noise-floor correction. Headline DEI = treated.EI_main - control.EI_main.
   EI_noise is reported as a within-sample diagnostic only.
 - QA-H9 fix: CpG sites already excluded from non-TpC set upstream (build_tpc_cds_bed.py).
 - Strand handling: '+' strand TpC -> ref=C, edited dir = C->T, noise = (A+G)/2 at C-ref.
   '-' strand TpC (i.e. GpA on + strand) -> ref=G, edited dir = G->A, noise = (C+T)/2 at G-ref.

Inputs: per-chrom bcftools TSVs at /mnt/data/gold_index/{BE4,Parent}_chr*.tsv
        cds_tpc.bed, cds_npc.bed (motif annotation)
Output: /mnt/data/dna_features/index/per_gene_dei.parquet
"""
import sys, os, pandas as pd, numpy as np
from collections import defaultdict

INDEX_DIR = '/mnt/data/gold_index'
TPC_BED = '/mnt/data/dna_features/index/cds_tpc.bed'
NPC_BED = '/mnt/data/dna_features/index/cds_npc.bed'
OUT_DIR = '/mnt/data/dna_features/index'
VAF_CAP = 0.40
MIN_POS_PER_GENE = 20

# --- Load motif annotation: (chrom, pos1based) -> (gene, motif_class, strand) ---
ann = {}
for bedp, mot in [(TPC_BED, 'tpc'), (NPC_BED, 'npc')]:
    n = 0
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom, s, e, gene, sc, strand = f
            ann[(chrom, int(s)+1)] = (gene, mot, strand)  # mpileup is 1-based
            n += 1
    print(f"loaded {n:,} {mot} positions from {bedp}", file=sys.stderr)
print(f"total annotated positions: {len(ann):,}", file=sys.stderr)

def parse_tsv_to_per_pos(path, label):
    """Parse bcftools-AD output. Returns dict (chrom,pos) -> dict with ref, A,C,G,T, tot, vaf_main, strand, motif, gene."""
    out = {}
    rows = 0
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5: continue
            rows += 1
            chrom = f[0]; pos = int(f[1]); ref = f[2].upper()
            key = (chrom, pos)
            if key not in ann: continue
            gene, motif, strand = ann[key]
            try:
                alts = f[3].split(',')
                counts = list(map(int, f[4].split(',')))
            except ValueError:
                continue
            if len(counts) != len(alts) + 1: continue
            # counts[0] = REF count; counts[1..] map to alts
            n_ref = counts[0]
            # tot excludes the <*> placeholder count (indels / ambiguous)
            tot = sum(counts) - (counts[-1] if alts[-1] == '<*>' else 0)
            if tot <= 0: continue
            ac = {'A':0,'C':0,'G':0,'T':0}
            ac[ref] = n_ref
            for a, c in zip(alts, counts[1:]):
                if a in ac: ac[a] = c
            # strand-specific main/noise
            if strand == '+':
                if ref != 'C': continue
                main = ac['T']
                noise = (ac['A'] + ac['G']) / 2.0
            else:
                if ref != 'G': continue
                main = ac['A']
                noise = (ac['C'] + ac['T']) / 2.0
            v_main = main / tot
            out[key] = dict(ref=ref, tot=tot, main=main, noise=noise, v_main=v_main,
                            strand=strand, motif=motif, gene=gene)
    print(f"  {label}: rows {rows:,} -> kept {len(out):,} annotated positions", file=sys.stderr)
    return out

CH = [f'chr{i}' for i in range(1,23)] + ['chrX']
all_T = {}; all_C = {}
for c in CH:
    be4p = f'{INDEX_DIR}/BE4_{c}.tsv'
    parp = f'{INDEX_DIR}/Parent_{c}.tsv'
    if not (os.path.exists(be4p) and os.path.exists(parp)):
        print(f"  skip {c}: missing", file=sys.stderr); continue
    print(f"=== {c} ===", file=sys.stderr)
    all_T.update(parse_tsv_to_per_pos(be4p, f'BE4 {c}'))
    all_C.update(parse_tsv_to_per_pos(parp, f'Parent {c}'))

# QA-M2 fix: ASYMMETRIC VAF filter. The previous "drop if VAF>cap in EITHER sample" killed
# real editing at high VAF in BE4 (clonal-expansion can fix editing events at 30-50% VAF in
# the clone). Now: drop a position iff the CONTROL (Parent) shows any meaningful variant,
# which indicates germline or pre-existing somatic — not BE-induced. Allow BE4 high VAF.
# QA-H5 fix: len(all_T) | len(all_C) was bitwise-OR (silent bug); use set union.
common = set(all_T) & set(all_C)
dropped_unpaired = len(set(all_T) | set(all_C)) - len(common)
dropped_germline = 0; whitelist = set()
PARENT_VAF_CAP = 0.05  # Parent-side: any detectable variant => not BE-specific
for key in common:
    if all_C[key]['v_main'] > PARENT_VAF_CAP:
        dropped_germline += 1
        continue
    whitelist.add(key)
print(f"\nposition filter: paired={len(common):,} | dropped (Parent VAF > {PARENT_VAF_CAP})={dropped_germline:,} | whitelist={len(whitelist):,}", file=sys.stderr)
print(f"unpaired positions dropped: {dropped_unpaired:,}", file=sys.stderr)

# Aggregate per gene (separately for treated and control, separately for TpC vs nonTpC)
agg = defaultdict(lambda: {  # gene -> dict
    'tT_pos':0, 'tT_tot':0, 'tT_main':0, 'tT_noise':0,
    'tC_pos':0, 'tC_tot':0, 'tC_main':0, 'tC_noise':0,
    'nT_pos':0, 'nT_tot':0, 'nT_main':0, 'nT_noise':0,
    'nC_pos':0, 'nC_tot':0, 'nC_main':0, 'nC_noise':0,
})
for key in whitelist:
    rec_T = all_T[key]; rec_C = all_C[key]
    g = rec_T['gene']; mot = rec_T['motif']
    p = 't' if mot == 'tpc' else 'n'
    d = agg[g]
    d[f'{p}T_pos'] += 1; d[f'{p}T_tot'] += rec_T['tot']; d[f'{p}T_main'] += rec_T['main']; d[f'{p}T_noise'] += rec_T['noise']
    d[f'{p}C_pos'] += 1; d[f'{p}C_tot'] += rec_C['tot']; d[f'{p}C_main'] += rec_C['main']; d[f'{p}C_noise'] += rec_C['noise']

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

# QA-H2 fix: NO double correction. Headline DEI = treated_main - control_main (rate-on-rate).
df['DEI_tpc'] = df.tT_EI_main - df.tC_EI_main
df['DEI_npc'] = df.nT_EI_main - df.nC_EI_main
# Within-sample noise sanity (treated only): how much higher is C->T than the A,G transversion floor?
df['noise_gap_tpc'] = df.tT_EI_main - df.tT_EI_noise
df['noise_gap_npc'] = df.nT_EI_main - df.nT_EI_noise

# coverage filter
df_full = df.copy()
df = df[(df.tT_n_pos >= MIN_POS_PER_GENE) & (df.tC_n_pos >= MIN_POS_PER_GENE)].copy()
df.to_parquet(f'{OUT_DIR}/per_gene_dei.parquet')

# ---- QC summary ----
def summ(a): return f'n={len(a)} mean={a.mean():.5g} med={a.median():.5g} q90={a.quantile(0.9):.5g}'
qc = [
    f"=== DNA Editing Index QC (post-fix) ===",
    f"VAF_CAP={VAF_CAP}  MIN_POS_PER_GENE={MIN_POS_PER_GENE}",
    f"shared paired positions: {len(common):,}",
    f"  dropped (high VAF in either): {dropped_high_vaf:,}  ({100*dropped_high_vaf/max(len(common),1):.1f}%)",
    f"  surviving whitelist:        {len(whitelist):,}",
    f"",
    f"genes after coverage filter: {len(df):,} of {len(df_full):,}",
    f"",
    f"--- Per-gene rates (low-VAF only) ---",
    f"Treated EI_main TpC:  {summ(df.tT_EI_main)}",
    f"Control EI_main TpC:  {summ(df.tC_EI_main)}",
    f"Treated EI_main nonTpC: {summ(df.nT_EI_main)}",
    f"Control EI_main nonTpC: {summ(df.nC_EI_main)}",
    f"",
    f"--- DNA EDITING INDEX (treated - control, rate) ---",
    f"DEI_tpc (HEADLINE): {summ(df.DEI_tpc)}",
    f"DEI_npc (motif-neg control, CpG-excluded): {summ(df.DEI_npc)}",
    f"",
    f"--- Sanity gates ---",
    f"[a] DEI_tpc mean > 0:                                {df.DEI_tpc.mean() > 0}  (need True)",
    f"[b] Mean(DEI_tpc) / Mean(DEI_npc):                    {df.DEI_tpc.mean() / max(df.DEI_npc.mean(), 1e-12):.2f}  (need >> 1)",
    f"[c] Fraction of genes with DEI_tpc > DEI_npc:        {(df.DEI_tpc > df.DEI_npc).mean():.3f}  (need > 0.70)",
    f"[d] Treated EI_main at TpC > Control EI_main at TpC: {(df.tT_EI_main > df.tC_EI_main).mean():.3f}",
    f"",
    f"--- Within-sample noise diagnostics (treated) ---",
    f"noise_gap at TpC (main - noise floor): {summ(df.noise_gap_tpc)}",
    f"noise_gap at nonTpC:                   {summ(df.noise_gap_npc)}",
]
qc_text = '\n'.join(qc)
print(qc_text)
with open(f'{OUT_DIR}/qc_summary.txt', 'w') as fh:
    fh.write(qc_text + '\n')
print(f"\nwrote {OUT_DIR}/per_gene_dei.parquet and qc_summary.txt", file=sys.stderr)
