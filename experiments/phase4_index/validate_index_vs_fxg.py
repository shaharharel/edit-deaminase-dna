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
    # Pearson of ranks = proper tie-corrected Spearman.
    # The shortcut formula 1 - 6 sum(d^2)/(n(n^2-1)) is WRONG when ties exist
    # (and our DEI has many ties at zero).
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else 0.0

print(f"genes in common: {len(m):,}")
# IMPORTANT: empirical DEI is a RATE per gene; compare to predicted MEAN (also a rate),
# NOT the sum which favors large genes. Use fxg_score_mean.
PRED_RATE = m.fxg_score_mean if 'fxg_score_mean' in m else m.fxg_score
PRED_SUM  = m.fxg_score if 'fxg_score' in m else PRED_RATE

print("\n=== Rate-vs-rate (the honest comparison) ===")
print(f"Spearman(empirical_DEI_tpc, predicted_fxg_MEAN):     {spearman(m.DEI_tpc, PRED_RATE):.3f}")
print(f"Spearman(empirical_DEI_npc, predicted_fxg_MEAN):     {spearman(m.DEI_npc, PRED_RATE):.3f}  <- motif-negative control")
print(f"Spearman(empirical_DEI_tpc, n_tpc_in_gene):          {spearman(m.DEI_tpc, m.n_tpc if 'n_tpc' in m else m.tT_n_pos):.3f}  <- gene-size / TpC-opportunity confound control")
print(f"Spearman(empirical_DEI_tpc, predicted_f_only):       {spearman(m.DEI_tpc, m.f_only):.3f}")
print(f"Spearman(empirical_DEI_tpc, predicted_g_only):       {spearman(m.DEI_tpc, m.g_only):.3f}")

print("\n=== Sum-vs-sum (gene-size driven; lower bar) ===")
print(f"Spearman(empirical_DEI_tpc * n_tpc, fxg_SUM):         {spearman(m.DEI_tpc * m.n_tpc, PRED_SUM) if 'n_tpc' in m else float('nan'):.3f}")

# top-K recall of high-DEI genes — using the rate
def recall_top_k(emp, pred, K=0.10):
    n = len(emp)
    top_emp = set(np.argsort(-emp.values)[:int(K*n)])
    top_pred = set(np.argsort(-pred.values)[:int(K*n)])
    return len(top_emp & top_pred) / len(top_emp)

print(f"\n=== Recall@top-K (rate-based ranking) ===")
print(f"Top-10% empirical DEI captured by top-10% predicted (rate): {recall_top_k(m.DEI_tpc, PRED_RATE):.3f}")
print(f"Top-5%  empirical DEI captured by top-5%  predicted (rate): {recall_top_k(m.DEI_tpc, PRED_RATE, 0.05):.3f}")
print(f"Top-1%  empirical DEI captured by top-1%  predicted (rate): {recall_top_k(m.DEI_tpc, PRED_RATE, 0.01):.3f}")

# enrichment in the motif-negative control (should be ~random)
print(f"\n=== Sanity: motif-negative control (should be ~K%) ===")
print(f"Top-10% empirical DEI_npc captured by top-10% predicted:    {recall_top_k(m.DEI_npc, PRED_RATE):.3f}")

m.to_parquet(OUT)
print(f"\nWrote {OUT}")
