#!/usr/bin/env python
"""Both Lj-BE clones failed the noise-floor criterion. The obvious next question is whether the
eA3A clones -- whose private-fraction 1.371x I have been quoting all night as a reference --
pass it. If they do not, that reference value is measuring the same thing.

Same criterion as the calibrator gate: alt>=1 & cov>=8 site count, and the fraction of those
below VAF 0.05. Node A's sound samples sit at 780k-1.14M and 51.8-67.2%.
"""
import numpy as np, glob, os

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
have = sorted({os.path.basename(p).replace("counts_", "").replace("_chr1.npz", "")
               for p in glob.glob(f"{FEAT}/counts_*_chr1.npz")})
print(f"  {'sample':26s} {'alt>=1':>12s} {'%<VAF.05':>9s} {'mean cov':>9s}  verdict")
for s in have:
    if len(glob.glob(f"{FEAT}/counts_{s}_chr*.npz")) < 23:
        print(f"  {s:26s} incomplete, skipped"); continue
    n = lo = ct = cn = 0
    for c in CH:
        D = np.load(f"{FEAT}/counts_{s}_chr{c}.npz")
        cov, alt = D["cov"], D["alt"]
        m = (alt >= 1) & (cov >= 8)
        n += int(m.sum()); ct += int(cov.sum()); cn += len(cov)
        v = alt[m] / np.maximum(cov[m], 1)
        lo += int((v < 0.05).sum())
    f = lo / max(n, 1)
    bad = []
    if f > 0.80: bad.append("sub-VAF")
    if n > 3 * 1_000_000: bad.append("count")     # 3x the ~1M healthy median on node A
    v = "FAILS " + "+".join(bad) if bad else "ok"
    print(f"  {s:26s} {n:>12,} {100*f:>8.1f}% {ct/max(cn,1):>9.2f}  {v}")
print("\n  healthy reference (node A): 780,594-1,135,489 sites, 51.8-67.2% sub-0.05")
print("  D10A-clone1, disqualified:  6,457,037 sites, 87.2% sub-0.05")
