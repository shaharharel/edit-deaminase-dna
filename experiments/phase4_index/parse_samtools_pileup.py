"""Parse samtools mpileup output (faster than bcftools mpileup-VCF pipeline).
Format: chrom pos ref depth bases quals
Computes per-gene DNA Editing Index with QA-fixed controls:
 - Asymmetric VAF cap (Parent only, 0.05) — preserve high-VAF BE4 editing
 - Shared whitelist: drop position from both samples if Parent shows variant
 - No double noise correction (DEI = treated_main - control_main)
 - CpG-excluded non-TpC motif-negative control
"""
import sys, os, pandas as pd, numpy as np
from collections import defaultdict
import re

OUTDIR = '/mnt/data/dna_features/index'
GOLD = '/mnt/data/gold_index'
TPC_BED = f'{OUTDIR}/cds_tpc.bed'
NPC_BED = f'{OUTDIR}/cds_npc.bed'
PARENT_VAF_CAP = 0.05
MIN_POS_PER_GENE = 20

# Load annotation: (chrom, pos1based) -> (gene, motif, strand)
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
    """Parse samtools mpileup bases string -> A/C/G/T counts."""
    out = {'A':0,'C':0,'G':0,'T':0}
    R = ref.upper()
    i = 0; n = len(bases)
    while i < n:
        c = bases[i]
        if c == '^':  # start-of-read followed by mapq char
            i += 2; continue
        if c == '$':
            i += 1; continue
        if c in '+-':  # indel: +N<seq> or -N<seq>
            j = i + 1
            while j < n and bases[j].isdigit(): j += 1
            try:
                ln_indel = int(bases[i+1:j])
                i = j + ln_indel
            except ValueError:
                i += 1
            continue
        if c == '*':
            i += 1; continue
        if c in '.,':
            if R in out: out[R] += 1
        elif c.upper() in 'ACGT':
            out[c.upper()] += 1
        i += 1
    return out

def parse_pileup_to_dict(path, label):
    """Return dict (chrom, pos) -> {ref, tot, main, noise, v_main, strand, motif, gene}"""
    out = {}
    n_rows = 0; n_skip_unann = 0; n_skip_ref = 0
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5: continue
            n_rows += 1
            chrom = f[0]; pos = int(f[1]); ref = f[2].upper()
            depth = int(f[3])
            if depth < 1: continue
            bases = f[4] if len(f) > 4 else ''
            key = (chrom, pos)
            if key not in ann: n_skip_unann += 1; continue
            gene, motif, strand = ann[key]
            ac = parse_bases(bases, ref)
            tot = ac['A'] + ac['C'] + ac['G'] + ac['T']
            if tot < 1: continue
            if strand == '+':
                if ref != 'C':
                    n_skip_ref += 1; continue
                main = ac['T']
                noise = (ac['A'] + ac['G']) / 2.0
            else:
                if ref != 'G':
                    n_skip_ref += 1; continue
                main = ac['A']
                noise = (ac['C'] + ac['T']) / 2.0
            out[key] = dict(tot=tot, main=main, noise=noise, v_main=main/tot,
                            strand=strand, motif=motif, gene=gene)
            if n_rows % 1_000_000 == 0:
                print(f"  {label}: parsed {n_rows:,} rows -> {len(out):,} kept", file=sys.stderr)
    print(f"  {label}: total {n_rows:,} rows | kept {len(out):,} | skipped unannot {n_skip_unann:,} | ref-mismatch {n_skip_ref:,}", file=sys.stderr)
    return out

print("Parsing BE4...", file=sys.stderr)
T = parse_pileup_to_dict(f'{GOLD}/BE4_samtools.mpileup', 'BE4')
print("Parsing Parent...", file=sys.stderr)
C = parse_pileup_to_dict(f'{GOLD}/Parent_samtools.mpileup', 'Parent')

# Asymmetric VAF whitelist: drop position iff Parent VAF > cap (control has variant => not BE-induced)
common = set(T) & set(C)
unpaired = len(set(T) | set(C)) - len(common)
dropped_germline = 0; whitelist = set()
for key in common:
    if C[key]['v_main'] > PARENT_VAF_CAP:
        dropped_germline += 1
        continue
    whitelist.add(key)
print(f"\nposition filter: paired={len(common):,} | dropped (Parent VAF > {PARENT_VAF_CAP})={dropped_germline:,} | whitelist={len(whitelist):,}", file=sys.stderr)
print(f"unpaired positions dropped: {unpaired:,}", file=sys.stderr)

# Per-gene aggregation
agg = defaultdict(lambda: {
    'tT_pos':0,'tT_tot':0,'tT_main':0,'tT_noise':0,
    'tC_pos':0,'tC_tot':0,'tC_main':0,'tC_noise':0,
    'nT_pos':0,'nT_tot':0,'nT_main':0,'nT_noise':0,
    'nC_pos':0,'nC_tot':0,'nC_main':0,'nC_noise':0,
})
for key in whitelist:
    rT = T[key]; rC = C[key]
    g = rT['gene']; mot = rT['motif']
    p = 't' if mot == 'tpc' else 'n'
    d = agg[g]
    d[f'{p}T_pos'] += 1; d[f'{p}T_tot'] += rT['tot']; d[f'{p}T_main'] += rT['main']; d[f'{p}T_noise'] += rT['noise']
    d[f'{p}C_pos'] += 1; d[f'{p}C_tot'] += rC['tot']; d[f'{p}C_main'] += rC['main']; d[f'{p}C_noise'] += rC['noise']

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
df.to_parquet(f'{OUTDIR}/per_gene_dei.parquet')

# QC summary
def s(a): return f'n={len(a)} mean={a.mean():.5g} med={a.median():.5g} q90={a.quantile(0.9):.5g}'
qc = [
    f"=== DNA Editing Index QC (samtools mpileup, post-fix) ===",
    f"PARENT_VAF_CAP={PARENT_VAF_CAP}  MIN_POS_PER_GENE={MIN_POS_PER_GENE}",
    f"paired positions: {len(common):,}",
    f"  dropped (germline / Parent variant): {dropped_germline:,} ({100*dropped_germline/max(len(common),1):.2f}%)",
    f"  surviving whitelist: {len(whitelist):,}",
    f"",
    f"genes after coverage filter: {len(df):,} of {len(df_full):,}",
    f"",
    f"--- Per-gene rates (low-VAF only) ---",
    f"Treated EI_main TpC:    {s(df.tT_EI_main)}",
    f"Control EI_main TpC:    {s(df.tC_EI_main)}",
    f"Treated EI_main nonTpC: {s(df.nT_EI_main)}",
    f"Control EI_main nonTpC: {s(df.nC_EI_main)}",
    f"",
    f"--- DNA EDITING INDEX (treated - control, rate) ---",
    f"DEI_tpc (HEADLINE): {s(df.DEI_tpc)}",
    f"DEI_npc (motif-NEG control, CpG-excluded): {s(df.DEI_npc)}",
    f"",
    f"--- Pre-registered sanity gates ---",
    f"[a] DEI_tpc mean > 0:                                {df.DEI_tpc.mean() > 0}",
    f"[b] Mean(DEI_tpc) / Mean(DEI_npc):                    {df.DEI_tpc.mean() / max(df.DEI_npc.mean(), 1e-12):.2f}  (need > 2x = convincing; > 1.3x = weak)",
    f"[c] Fraction genes DEI_tpc > DEI_npc:                {(df.DEI_tpc > df.DEI_npc).mean():.3f}  (need > 0.70)",
    f"[d] Fraction genes Treated TpC > Control TpC:        {(df.tT_EI_main > df.tC_EI_main).mean():.3f}",
    f"",
    f"--- Within-sample noise diagnostics (treated) ---",
    f"noise_gap at TpC (main - noise floor): {s(df.noise_gap_tpc)}",
    f"noise_gap at nonTpC:                   {s(df.noise_gap_npc)}",
]
text = '\n'.join(qc)
print(text)
with open(f'{OUTDIR}/qc_summary.txt', 'w') as fh: fh.write(text + '\n')
print(f"\nwrote {OUTDIR}/per_gene_dei.parquet and qc_summary.txt", file=sys.stderr)
