"""Per-gene f x g prediction from rAPOBEC1 motif spectrum + 1Mb-bin accessibility.
For each TpC site in cds_tpc.bed:
  f_site = trinuc log-odds from rAPOBEC1 sibling sources (computed against random-C bg)
  g_site = z-scored 1Mb-bin accessibility composite (mean of DNase / ATAC / R-loop / H3K27ac, z-scored)
  score_site = z(f_site) + z(g_site)
Per-gene aggregate (3 versions):
  - mean of site scores
  - sum of site scores (rate per gene proxy)
  - top-K-mean (concentration at the strongest sites)
Outputs per_gene_fxg.parquet  with cols: gene, n_tpc, f_only, g_only, fxg_score, fxg_score_mean
"""
import sys, pandas as pd, numpy as np
from pyfaidx import Fasta
from collections import Counter, defaultdict

CDS_TPC_BED = sys.argv[1]
HG38 = sys.argv[2]
BINS_PARQUET = sys.argv[3]    # bins_1mb_v3.parquet  (with chrom,bin,dnase,atac,rloop,h3k27ac)
CANONICAL = sys.argv[4]       # canonical_sites.parquet (for rAPOBEC1 spectrum estimation)
OUT = sys.argv[5]

fa = Fasta(HG38)
W = 10
COMP = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))

# --- 1) rAPOBEC1 trinucleotide log-odds spectrum from sibling sources ---
canon = pd.read_parquet(CANONICAL)
rapob = canon[canon.family == 'rAPOBEC1']
pos_tn = Counter()
for r in rapob.itertuples():
    c = r.context
    if len(c) == 2*W+1 and c[W] == 'C':
        pos_tn[c[W-1] + c[W] + c[W+1]] += 1
# background trinuc freq: sample 200k random genomic Cs
bg_tn = Counter()
CH = [f'chr{i}' for i in range(1,23)] + ['chrX']
import numpy.random as nr
rng = nr.default_rng(0)
CL = {c: len(fa[c]) for c in CH}
wt = np.array([CL[c] for c in CH], float); wt /= wt.sum()
n = 0
while n < 200_000:
    c = rng.choice(CH, p=wt); p = int(rng.integers(W+2, CL[c]-W-2))
    try: w = fa[c][p-1-W:p+W].seq.upper()
    except: continue
    if len(w) != 2*W+1 or 'N' in w: continue
    if w[W] == 'G': w = rc(w)
    if w[W] != 'C': continue
    bg_tn[w[W-1]+w[W]+w[W+1]] += 1; n += 1
tot_p = sum(pos_tn.values()); tot_b = sum(bg_tn.values())
spec = {t: np.log((pos_tn[t]+1)/(tot_p+1)) - np.log((bg_tn[t]+1)/(tot_b+1))
        for t in set(list(pos_tn)+list(bg_tn))}
print(f"rAPOBEC1 spectrum estimated from {tot_p:,} positives. Top trinucs by log-odds:", file=sys.stderr)
for t in sorted(spec, key=lambda x: -spec[x])[:5]:
    print(f"   {t}: {spec[t]:+.3f}", file=sys.stderr)

# --- 2) accessibility composite per 1Mb bin ---
b = pd.read_parquet(BINS_PARQUET)
for col in ['dnase','atac','rloop','h3k27ac']:
    if col not in b.columns:
        raise SystemExit(f"bins file missing column {col}")
g_arr = b[['dnase','atac','rloop','h3k27ac']].values
mu = np.nanmean(g_arr, 0); sd = np.nanstd(g_arr, 0) + 1e-6
g_arr_z = (g_arr - mu) / sd
g_score = g_arr_z.mean(axis=1)  # composite g per bin
acc = {(r.chrom, int(r.bin)): g_score[i] for i,r in enumerate(b.itertuples())}

# --- 3) score every CDS TpC site, aggregate per gene ---
per_gene_f = defaultdict(list)
per_gene_g = defaultdict(list)
per_gene_fxg = defaultdict(list)
seen = 0
with open(CDS_TPC_BED) as fh:
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 6: continue
        chrom, s, e, gene, _, strand = f
        pos = int(s)  # 0-based; the C itself
        try:
            w = fa[chrom][pos-W:pos+W+1].seq.upper()
        except: continue
        if len(w) != 2*W+1 or 'N' in w: continue
        if strand == '-':
            w = rc(w)
        if w[W] != 'C': continue
        tn = w[W-1] + w[W] + w[W+1]
        f_site = spec.get(tn, -5.0)
        g_site = acc.get((chrom, pos // 1_000_000), 0.0)
        per_gene_f[gene].append(f_site)
        per_gene_g[gene].append(g_site)
        seen += 1
print(f"scored {seen:,} CDS TpC sites across {len(per_gene_f):,} genes", file=sys.stderr)

# z-score the per-site scores across the WHOLE distribution, then aggregate per gene
all_f = np.concatenate([np.array(v) for v in per_gene_f.values()])
all_g = np.concatenate([np.array(v) for v in per_gene_g.values()])
fm, fs = all_f.mean(), all_f.std() + 1e-6
gm, gs = all_g.mean(), all_g.std() + 1e-6

rows = []
for gene in per_gene_f:
    fv = (np.array(per_gene_f[gene]) - fm) / fs
    gv = (np.array(per_gene_g[gene]) - gm) / gs
    fxg = fv + gv
    rows.append({
        'gene': gene,
        'n_tpc': len(fv),
        'f_only': float(fv.mean()),
        'g_only': float(gv.mean()),
        'fxg_score_mean': float(fxg.mean()),
        'fxg_score': float(fxg.sum()),   # sum ~ per-gene predicted burden
        'fxg_top10pct_mean': float(np.sort(fxg)[-max(1,len(fxg)//10):].mean()),
    })
out = pd.DataFrame(rows)
out.to_parquet(OUT)
print(f"wrote {OUT}  ({len(out):,} genes)", file=sys.stderr)
