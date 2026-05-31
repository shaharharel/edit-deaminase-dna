"""Alternative DEI: count discrete editing-candidate POSITIONS per gene
(positions where variant VAF is in [VAF_MIN, VAF_MAX] = above sequencing noise but below germline).
This matches the user's intuition: 1 positive read at one position counts as 1 'edit event' for the gene.
Compare counts at TpC vs non-TpC, and treated vs control.
"""
import sys, os, pandas as pd, numpy as np
from collections import defaultdict

INDEX_DIR = '/mnt/data/gold_index'
TPC_BED = '/mnt/data/dna_features/index/cds_tpc.bed'
NPC_BED = '/mnt/data/dna_features/index/cds_npc.bed'
OUT_DIR = '/mnt/data/dna_features/index'
VAF_MIN = 0.02   # above sequencing error rate (~1%)
VAF_MAX = 0.30   # below germline het (~50%); 0.30 is conservative
MIN_POS_PER_GENE = 20

ann = {}
for bedp, mot in [(TPC_BED, 'tpc'), (NPC_BED, 'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom, s, e, gene, sc, strand = f
            ann[(chrom, int(s)+1)] = (gene, mot, strand)
print(f"loaded {len(ann):,} annotated positions", file=sys.stderr)

def parse_to_dict(path, label):
    """Per position: (vaf_main, gene, motif, strand, depth_ok)."""
    out = {}
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5: continue
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
            n_ref = counts[0]
            tot = sum(counts) - (counts[-1] if alts[-1] == '<*>' else 0)
            if tot < 10: continue   # require min depth for VAF stability
            ac = {'A':0,'C':0,'G':0,'T':0}
            ac[ref] = n_ref
            for a, c in zip(alts, counts[1:]):
                if a in ac: ac[a] = c
            if strand == '+':
                if ref != 'C': continue
                main = ac['T']
            else:
                if ref != 'G': continue
                main = ac['A']
            out[key] = (main / tot, gene, motif, strand)
    return out

CH = [f'chr{i}' for i in range(1,23)] + ['chrX']
T = {}; C = {}
for c in CH:
    p1 = f'{INDEX_DIR}/BE4_{c}.tsv'
    p2 = f'{INDEX_DIR}/Parent_{c}.tsv'
    if not (os.path.exists(p1) and os.path.exists(p2)): continue
    print(f"{c}: parsing", file=sys.stderr)
    T.update(parse_to_dict(p1, 'BE4'))
    C.update(parse_to_dict(p2, 'Parent'))
print(f"BE4 positions: {len(T):,}  Parent positions: {len(C):,}", file=sys.stderr)

# per-gene counts: edited positions in band [VAF_MIN, VAF_MAX], at TpC vs nonTpC
agg = defaultdict(lambda: {'tT_e':0,'tT_n':0,'tC_e':0,'tC_n':0,
                            'nT_e':0,'nT_n':0,'nC_e':0,'nC_n':0})
both = set(T) & set(C)
for k in both:
    vt, gene, motif, _ = T[k]
    vc = C[k][0]
    p = 't' if motif == 'tpc' else 'n'
    d = agg[gene]
    d[f'{p}T_n'] += 1
    d[f'{p}C_n'] += 1
    if VAF_MIN <= vt <= VAF_MAX: d[f'{p}T_e'] += 1
    if VAF_MIN <= vc <= VAF_MAX: d[f'{p}C_e'] += 1

rows = []
for gene, d in agg.items():
    row = {'gene': gene}
    for lbl in ('tT','tC','nT','nC'):
        row[f'{lbl}_n'] = d[f'{lbl}_n']
        row[f'{lbl}_e'] = d[f'{lbl}_e']
        row[f'{lbl}_rate'] = d[f'{lbl}_e'] / max(d[f'{lbl}_n'], 1)
    rows.append(row)
df = pd.DataFrame(rows)
df['DEI_tpc'] = df.tT_rate - df.tC_rate   # treated minus control rate (positions per total)
df['DEI_npc'] = df.nT_rate - df.nC_rate
df_full = df.copy()
df = df[(df.tT_n >= MIN_POS_PER_GENE) & (df.tC_n >= MIN_POS_PER_GENE)].copy()
df.to_parquet(f'{OUT_DIR}/per_gene_dei_count.parquet')

def s(a): return f'n={len(a)} mean={a.mean():.4g} med={a.median():.4g} q90={a.quantile(0.9):.4g}'
print(f"=== Position-count DEI (VAF in [{VAF_MIN},{VAF_MAX}]) ===")
print(f"surviving genes: {len(df):,} of {len(df_full):,}")
print(f"treated TpC edit-pos rate: {s(df.tT_rate)}")
print(f"control TpC edit-pos rate: {s(df.tC_rate)}")
print(f"treated nonTpC edit-pos rate: {s(df.nT_rate)}")
print(f"control nonTpC edit-pos rate: {s(df.nC_rate)}")
print()
print(f"DEI_tpc (treated-control at TpC): {s(df.DEI_tpc)}")
print(f"DEI_npc (treated-control at nonTpC): {s(df.DEI_npc)}")
print()
print(f"Mean(DEI_tpc)/Mean(DEI_npc): {df.DEI_tpc.mean()/max(df.DEI_npc.mean(),1e-12):.2f}")
print(f"Fraction DEI_tpc > DEI_npc: {(df.DEI_tpc > df.DEI_npc).mean():.3f}")
print(f"Treated TpC count > Treated nonTpC count: {(df.tT_e > df.nT_e).mean():.3f}")
print()
print(f"--- Treated EVENT counts per gene ---")
print(f"  TpC edit events: total {df.tT_e.sum():,}; median {df.tT_e.median():.0f}; q90 {df.tT_e.quantile(0.9):.0f}")
print(f"  nonTpC edit events: total {df.nT_e.sum():,}; median {df.nT_e.median():.0f}; q90 {df.nT_e.quantile(0.9):.0f}")
print(f"  Total ratio TpC/nonTpC events: {df.tT_e.sum()/max(df.nT_e.sum(),1):.2f}")
print(f"  Per-position rate: TpC {df.tT_e.sum()/df.tT_n.sum():.5g}  nonTpC {df.nT_e.sum()/df.nT_n.sum():.5g}")
print(f"  Ratio of per-position rates: {(df.tT_e.sum()/df.tT_n.sum())/max(df.nT_e.sum()/df.nT_n.sum(),1e-12):.3f}")
