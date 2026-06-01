"""Validate empirical DNA Editing Index against predicted f x g per gene.

QA-driven design:
 - Empirical DEI is a RATE per gene; compare to predicted RATE (fxg_score_mean), not sum.
 - Spearman: Pearson-of-ranks (proper tie correction).
 - Pre-registered headline gate: motif-specificity ratio = Spearman(DEI_tpc, fxg) / Spearman(DEI_npc, fxg).
   Convincing: > 2x (with Spearman > 0.10). Weak: 1.3-2x. Honest negative: < 1.3x.
 - Partial Spearman controlling for n_tpc (gene-size confound).
 - Recall@K computed only over genes with DEI > 0 (zero-DEI ties make full-pool recall meaningless).
 - fxg_top10pct_mean as a secondary (rank by the most-likely-edited sites in each gene).
"""
import sys, pandas as pd, numpy as np
EMP = sys.argv[1]; PRED = sys.argv[2]; OUT = sys.argv[3]
e = pd.read_parquet(EMP); p = pd.read_parquet(PRED)
m = e.merge(p, on='gene', how='inner')
# drop NaN
m = m.dropna(subset=['DEI_tpc','DEI_npc','fxg_score_mean','f_only','g_only']).copy()
print(f"genes in common (after NaN drop): {len(m):,}")

def spearman(x, y):
    rx = pd.Series(x).rank().values; ry = pd.Series(y).rank().values
    rx = rx - rx.mean(); ry = ry - ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else 0.0

def partial_spearman(x, y, z):
    """Spearman(x, y | z) via Pearson on residuals of rank(x) and rank(y) regressed on rank(z)."""
    rx = pd.Series(x).rank().values; ry = pd.Series(y).rank().values; rz = pd.Series(z).rank().values
    def resid(a, b):
        b = b - b.mean(); a = a - a.mean()
        return a - ((a*b).sum() / max((b*b).sum(), 1e-12)) * b
    ex, ey = resid(rx, rz), resid(ry, rz)
    denom = np.sqrt((ex**2).sum() * (ey**2).sum())
    return float((ex*ey).sum() / denom) if denom > 0 else 0.0

print("\n=== HEADLINE (rate-vs-rate; motif specificity) ===")
sp_tpc = spearman(m.DEI_tpc, m.fxg_score_mean)
sp_npc = spearman(m.DEI_npc, m.fxg_score_mean)
sp_top = spearman(m.DEI_tpc, m.fxg_top10pct_mean) if 'fxg_top10pct_mean' in m else float('nan')
print(f"Spearman(DEI_tpc, fxg_mean):        {sp_tpc:+.4f}")
print(f"Spearman(DEI_npc, fxg_mean):        {sp_npc:+.4f}  <- motif-NEG control (should be near 0)")
print(f"Spearman(DEI_tpc, fxg_top10pct):    {sp_top:+.4f}  <- alt: max-site signal")
ratio = sp_tpc / sp_npc if abs(sp_npc) > 1e-3 else float('inf')
diff  = sp_tpc - sp_npc
print(f"MOTIF-SPECIFICITY RATIO (tpc/npc):  {ratio:+.2f}x")
print(f"MOTIF-SPECIFICITY DIFFERENCE:       {diff:+.4f}")

print("\n=== Decomposition ===")
print(f"Spearman(DEI_tpc, f_only):          {spearman(m.DEI_tpc, m.f_only):+.4f}")
print(f"Spearman(DEI_tpc, g_only):          {spearman(m.DEI_tpc, m.g_only):+.4f}")

print("\n=== Confound controls ===")
print(f"Spearman(DEI_tpc, n_tpc):           {spearman(m.DEI_tpc, m.n_tpc):+.4f}  <- gene-size confound (should be small)")
print(f"Spearman(fxg_mean, n_tpc):          {spearman(m.fxg_score_mean, m.n_tpc):+.4f}  <- predicted-side gene-size confound")
ps = partial_spearman(m.DEI_tpc, m.fxg_score_mean, m.n_tpc)
ps_npc = partial_spearman(m.DEI_npc, m.fxg_score_mean, m.n_tpc)
print(f"Partial Spearman(DEI_tpc, fxg_mean | n_tpc):  {ps:+.4f}  <- clean signal after partialling gene size")
print(f"Partial Spearman(DEI_npc, fxg_mean | n_tpc):  {ps_npc:+.4f}  <- motif-NEG with same control")

print("\n=== Recall (restricted to non-zero DEI; zero-DEI ties make full-pool meaningless) ===")
nz = m[m.DEI_tpc > 0].copy()
print(f"  non-zero DEI_tpc genes: {len(nz):,} / {len(m):,}")
def recall_at(emp, pred, K):
    if len(emp) == 0: return float('nan')
    o_emp = np.argsort(-emp.values); o_pred = np.argsort(-pred.values)
    k = max(1, int(K*len(emp)))
    top_e = set(o_emp[:k]); top_p = set(o_pred[:k])
    return len(top_e & top_p) / k
for K in (0.10, 0.05, 0.01):
    print(f"  R@{int(K*100)}% (DEI_tpc>0 vs fxg_mean):  {recall_at(nz.DEI_tpc, nz.fxg_score_mean, K):.3f}  (random ~ {K:.2f})")

print("\n=== Pre-registered gate ===")
print(f"Headline ratio (Spearman_tpc / Spearman_npc): {ratio:+.2f}x")
if abs(sp_tpc) > 0.10 and ratio > 2.0:
    verdict = "CONVINCING signal (DNA at single-clone 50x)"
elif abs(sp_tpc) > 0.05 and ratio > 1.3:
    verdict = "WEAK signal (suggestive, needs replication)"
else:
    verdict = "HONEST NEGATIVE (consistent with depth ceiling)"
print(f"Verdict: {verdict}")

m.to_parquet(OUT)
print(f"\nWrote {OUT}")
