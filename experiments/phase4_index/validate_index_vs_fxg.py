"""Validate empirical DNA Editing Index (Doman BE4-Parent) against the f x g model prediction.
Per gene:
  empirical  = DEI_tpc (treated - control, low-VAF, mismatch-direction corrected, at TpC)
  predicted  = aggregate f(motif spectrum) x g(accessibility) score over the gene's TpC sites
Spearman + per-tier (Tier A cancer drivers, etc.) breakdown.
"""
import sys, pandas as pd, numpy as np

EMP = sys.argv[1]   # per_gene_index.parquet from compute_dna_editing_index.py
PRED = sys.argv[2]  # per_gene_fxg.parquet from predict_per_gene_fxg.py (separate)
OUT = sys.argv[3]

e = pd.read_parquet(EMP)
p = pd.read_parquet(PRED)
m = e.merge(p, on='gene', how='inner')

def spearman(x, y):
    n = len(x)
    rx = pd.Series(x).rank(); ry = pd.Series(y).rank()
    d = rx - ry
    return 1 - 6 * (d**2).sum() / (n * (n**2 - 1))

print(f"genes in common: {len(m):,}")
print(f"Spearman(empirical_DEI_tpc, predicted_fxg):     {spearman(m.DEI_tpc, m.fxg_score):.3f}")
print(f"Spearman(empirical_DEI_npc, predicted_fxg):     {spearman(m.DEI_npc, m.fxg_score):.3f}  <- motif-negative control")
print(f"Spearman(empirical_DEI_tpc, gene_density):      {spearman(m.DEI_tpc, m.gene_density if 'gene_density' in m else m.tT_n_pos):.3f}  <- confound control")
print(f"Spearman(empirical_DEI_tpc, predicted_f_only):  {spearman(m.DEI_tpc, m.f_only):.3f}")
print(f"Spearman(empirical_DEI_tpc, predicted_g_only):  {spearman(m.DEI_tpc, m.g_only):.3f}")

# top-K recall of high-DEI genes
def recall_top_k(emp, pred, K=0.10):
    n = len(emp)
    top_emp = set(np.argsort(-emp.values)[:int(K*n)])
    top_pred = set(np.argsort(-pred.values)[:int(K*n)])
    return len(top_emp & top_pred) / len(top_emp)

print(f"\nRecall of top-10% empirically-edited genes within top-10% predicted: {recall_top_k(m.DEI_tpc, m.fxg_score):.3f}")
print(f"Recall of top-1% empirically-edited genes within top-1% predicted:   {recall_top_k(m.DEI_tpc, m.fxg_score, 0.01):.3f}")

m.to_parquet(OUT)
