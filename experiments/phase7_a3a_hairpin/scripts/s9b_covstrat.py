#!/usr/bin/env python
"""6(d) part 2. Part 1 found a REAL imbalance: the positive rate slides 1.322x -> 0.621x
across Parent-coverage deciles (a 2.1x gradient), and hairpin sites sit at slightly lower
coverage (mean stem 1.716 in the lowest decile vs 1.651 in the highest; frac cov<10 rises
0.029 -> 0.071 from stem>=0 to stem>=8). That is a live confounding path:
    hairpin -> lower mappability -> lower coverage -> higher positive rate
So the confound is not hypothetical and it had to be measured, not asserted away.

The decisive question is not whether the path exists but how much of the hairpin signal it
carries. Answer it by STRATIFICATION, which needs no retraining: if the hairpin positive-rate
ratio survives inside every coverage decile, the gradient cannot be its source. This is the
same instrument that cleared GC (10/10 deciles) and complexity.

Also reports the Mantel-Haenszel pooled estimate -- the stratified summary that a decile
table implies but does not state.
"""
import numpy as np

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
d = np.load(f"{FEAT}/a3a_trainset_v3.npz", allow_pickle=True)
chrom, pos, y, stem = d["chrom"], d["pos"], d["y"], d["stem"]

cov = np.full(len(y), -1, dtype=np.int32)
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    upos = U["pos1"] if "pos1" in U.files else U["pos"] + 1
    m = chrom == c
    if not m.any(): continue
    idx = np.clip(np.searchsorted(upos, pos[m]), 0, len(upos) - 1)
    ok = upos[idx] == pos[m]
    cov[np.flatnonzero(m)[ok]] = P["cov"][idx[ok]]
    del U, P
assert (cov >= 0).all(), "join incomplete"

qs = np.percentile(cov, np.arange(0, 101, 10))
for S_MIN in (6, 8):
    hp = stem >= S_MIN
    base_all = (y[hp] == 1).mean() / (y == 1).mean()
    print(f"\n=== stem>={S_MIN}: hairpin positive-rate ratio, STRATIFIED by Parent coverage ===")
    print(f"  unstratified: {base_all:.3f}x   (n_hp={int(hp.sum()):,})")
    print(f"  {'decile':>7s} {'cov':>10s} {'n_hp':>8s} {'p_hp':>8s} {'p_flat':>8s} {'ratio':>7s}")
    num = den = 0.0
    for i in range(10):
        lo, hi = qs[i], qs[i + 1]
        m = (cov >= lo) & (cov <= hi) if i == 9 else (cov >= lo) & (cov < hi)
        a = m & hp
        b = m & ~hp
        if a.sum() < 200: continue
        p_hp = (y[a] == 1).mean(); p_fl = (y[b] == 1).mean()
        # Mantel-Haenszel on the 2x2 (hairpin/flat) x (pos/neg) within this stratum
        A = float((y[a] == 1).sum()); B = float((y[a] == 0).sum())
        Cc = float((y[b] == 1).sum()); D = float((y[b] == 0).sum())
        T = A + B + Cc + D
        num += A * D / T; den += B * Cc / T
        print(f"  {i:>7d} {lo:4.0f}-{hi:<5.0f} {int(a.sum()):>8,} {p_hp:>8.5f} {p_fl:>8.5f} "
              f"{p_hp/p_fl:>6.3f}x")
    print(f"  MANTEL-HAENSZEL pooled OR (coverage-adjusted): {num/den:.3f}")
    print(f"  crude OR for comparison: "
          f"{((y[hp]==1).sum()/(y[hp]==0).sum())/((y[~hp]==1).sum()/(y[~hp]==0).sum()):.3f}")
