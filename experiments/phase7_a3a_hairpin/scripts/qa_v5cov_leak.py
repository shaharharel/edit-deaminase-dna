#!/usr/bin/env python
"""Check 1 on a3a_trainset_v5cov.npz -- the set that produced last tick's 5.119x.

I reported that number as putting the transfer test's positive control on controlled
ground. A result I have just leaned on is exactly the one to audit. Bug family: a
convention correct locally, consumed downstream as universal. v5cov was built by
RE-KEYING v5's negatives on coverage, so every join is a chance to reintroduce one.

Sampling note: this file is BLOCK-ORDERED (positives then negatives), so a head slice
would sample ~45% positives instead of 9.4%. Sample RANDOMLY and ASSERT the base rate.
"""
import numpy as np

P = "/data/a3a/feat/a3a_trainset_v5cov.npz"
d = np.load(P, allow_pickle=True)
y = d["y"].astype(bool)
win, tri, strand, cov = d["win"], d["tri"], d["strand"], d["cov"]
n = len(y)
print("=== base rate, ceiling, and the sampling assertion ===")
br = y.mean()
print(f"  n={n:,}  pos={int(y.sum()):,}  base_rate={br:.4f}  ceiling={1/br:.3f}x")
rng = np.random.default_rng(0)
idx = rng.choice(n, size=min(200_000, n), replace=False)
print(f"  random sample base rate {y[idx].mean():.4f} (must track {br:.4f}, not ~0.45)")
assert abs(y[idx].mean() - br) < 0.01, "SAMPLING BUG -- block-ordered slice"

W = win.shape[1] if win.ndim == 2 else None
print(f"\n=== window shape {win.shape}, dtype {win.dtype} ===")
c = W // 2 if W else None

print("\n=== CHECK 1a: focal base must be C for 100% of BOTH classes ===")
foc = win[:, c]
for lab, m in (("positives", y), ("negatives", ~y)):
    vals, cnts = np.unique(foc[m], return_counts=True)
    top = dict(zip(vals.tolist(), cnts.tolist()))
    frac = {int(k): v / m.sum() for k, v in top.items()}
    print(f"  {lab:10s} focal base distribution: " +
          " ".join(f"{k}:{v:.4f}" for k, v in sorted(frac.items())))

print("\n=== CHECK 1b: ACGT at -1 and +1 must be IDENTICAL between classes ===")
for off, lab in ((-1, "-1"), (1, "+1")):
    j = c + off
    rows = []
    for cl, m in (("pos", y), ("neg", ~y)):
        v, k = np.unique(win[m, j], return_counts=True)
        p = dict(zip(v.tolist(), (k / m.sum()).tolist()))
        rows.append((cl, p))
    keys = sorted(set(rows[0][1]) | set(rows[1][1]))
    dmax = max(abs(rows[0][1].get(k, 0) - rows[1][1].get(k, 0)) for k in keys)
    print(f"  offset {lab}:  " + "  ".join(
        f"{k}: pos {rows[0][1].get(k,0):.4f} / neg {rows[1][1].get(k,0):.4f}" for k in keys))
    print(f"           max |pos-neg| = {dmax:.4f}  {'OK' if dmax < 0.02 else '*** MISMATCH ***'}")

print("\n=== CHECK: trinuc and strand matching PRESERVED by the coverage re-key ===")
for name, arr in (("trinuc", tri), ("strand", strand)):
    vp, kp = np.unique(arr[y], return_counts=True)
    vn, kn = np.unique(arr[~y], return_counts=True)
    pp = dict(zip(vp.tolist(), (kp / y.sum()).tolist()))
    pn = dict(zip(vn.tolist(), (kn / (~y).sum()).tolist()))
    keys = sorted(set(pp) | set(pn))
    dmax = max(abs(pp.get(k, 0) - pn.get(k, 0)) for k in keys)
    print(f"  {name:7s} " + "  ".join(f"{k}: {pp.get(k,0):.4f}/{pn.get(k,0):.4f}" for k in keys))
    print(f"          max |pos-neg| = {dmax:.4f}  {'OK' if dmax < 0.02 else '*** NOT MATCHED ***'}")

print("\n=== CHECK: did coverage matching actually WORK? (that is the whole point) ===")
print(f"  cov present for {(cov>=0).mean():.4f} of rows")
qs = [10, 25, 50, 75, 90]
cp = np.percentile(cov[y & (cov >= 0)], qs)
cn = np.percentile(cov[~y & (cov >= 0)], qs)
print("  percentile   positives   negatives   diff")
for q, a, b in zip(qs, cp, cn):
    print(f"   {q:>3d}%       {a:8.2f}   {b:8.2f}   {a-b:+6.2f}")
mp, mn = cov[y & (cov >= 0)].mean(), cov[~y & (cov >= 0)].mean()
print(f"  mean         {mp:8.2f}   {mn:8.2f}   {mp-mn:+6.2f}  "
      f"({'MATCHED' if abs(mp-mn) < 1.0 else '*** STILL SKEWED ***'})")
