#!/usr/bin/env python
"""CHECK 6 applied to the headline number: is 5.036x carried by one chromosome?

The tail is pooled across five held-out folds. A pooled number can be produced by one fold
contributing almost all of it -- in which case the model is fitting something chromosome-
specific and the figure would not survive a different split. Nothing tonight has looked.

Reports the tail WITHIN each fold's own held-out set, with n for each, plus a leave-one-
fold-out recomputation of the pooled value. If dropping any single fold moves the pooled
number outside the ~12% empirical noise floor, the result depends on that fold.
"""
import numpy as np

FEAT = "/data/a3a/feat"
d = np.load(f"{FEAT}/a3a_trainset_v5cov.npz", allow_pickle=True)
y = d["y"].astype(float); chrom = d["chrom"]
p = np.load(f"{FEAT}/oof_v5cov_hairpin+sequence.npy")
if p.min() < 0 or p.max() > 1: p = 1 / (1 + np.exp(-p))
N = len(y); base = y.mean()
chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
folds = [chroms[i::5] for i in range(5)]
K = max(int(round(N * 0.001)), 1)     # derived from THIS file, never hardcoded: v5 is 924
print(f"n={N:,}  base rate {base:.5f}  top-0.1% is k={K}  pooled = "
      f"{y[np.argpartition(-p, K-1)[:K]].mean()/base:.3f}x")

rng = np.random.default_rng(0)
def tail(mask, pct=0.1):
    idx = np.flatnonzero(mask)
    k = max(int(round(len(idx) * pct / 100)), 1)
    sel = idx[np.argpartition(-p[idx], k - 1)[:k]]
    b = y[idx].mean()
    return y[sel].mean() / b, k, b

print(f"\n=== within each held-out fold (its own base rate, its own top 0.1%) ===")
print(f"  {'fold':>4} {'chroms':<22} {'n_test':>9} {'n_top':>6} {'base':>7} {'tail':>8}")
vals = []
for i, hold in enumerate(folds):
    m = np.isin(chrom, hold)
    e, k, b = tail(m)
    vals.append(e)
    print(f"  {i:>4} {','.join(hold):<22} {int(m.sum()):>9,} {k:>6,} {b:>7.4f} {e:>7.3f}x")
v = np.array(vals)
print(f"  spread {v.min():.3f}-{v.max():.3f}   mean {v.mean():.3f}   sd {v.std():.3f} "
      f"({100*v.std()/v.mean():.1f}% of mean)")

print(f"\n=== leave-one-fold-out: does the pooled number depend on any single fold? ===")
for i, hold in enumerate(folds):
    keep = ~np.isin(chrom, hold)
    e, k, b = tail(keep)
    print(f"  drop fold {i} ({','.join(hold):<20}) pooled {e:>6.3f}x  n_top={k:,}")
print(f"  noise floor is ~12%; a swing beyond that on one drop means the fold carries it.")

print(f"\n=== and per chromosome, the coarsest version of the same question ===")
rows = []
for c in chroms:
    m = chrom == c
    if m.sum() < 20000: continue
    e, k, b = tail(m)
    rows.append((e, c, int(m.sum()), k))
rows.sort()
print(f"  weakest 3: " + "  ".join(f"chr{c} {e:.2f}x (n={n:,})" for e, c, n, k in rows[:3]))
print(f"  strongest 3: " + "  ".join(f"chr{c} {e:.2f}x (n={n:,})" for e, c, n, k in rows[-3:]))
