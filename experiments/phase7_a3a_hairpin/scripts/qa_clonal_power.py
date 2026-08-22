#!/usr/bin/env python
"""CHECK 6 applied to the whole endpoint: can hairpin enrichment be measured on CLONAL
variants at all, with these data?

Last tick established that the specific-site sets are 91-94% below VAF 0.15, and that hairpin
enrichment reverses between VAF strata. The obvious next question -- does the enrichment
survive when restricted to genuinely clonal variants -- can only be answered if there ARE
enough clonal variants. Nobody has counted.

For each sample: specific sites at VAF>=0.35, and the EXPECTED number of hairpin sites among
them at stem>=8 (n x p_bg, p_bg ~ 0.002). If that expectation is a handful, the test is
powerless and saying so is the result.
"""
import sys, glob, os
import numpy as np

FEAT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
CUTS = [0.15, 0.25, 0.35]

samples = sorted({os.path.basename(p).replace("counts_", "").replace("_chr1.npz", "")
                  for p in glob.glob(f"{FEAT}/counts_*_chr1.npz")})
samples = [s for s in samples if s != "Parent" and
           len(glob.glob(f"{FEAT}/counts_{s}_chr*.npz")) >= 23]

# background hairpin rate on the Parent-silent set
bh = bn = 0
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15)
    bh += int((sil & (U["stem"] >= 8)).sum()); bn += int(sil.sum())
p_bg = bh / max(bn, 1)
print(f"  background P(stem>=8) on the Parent-silent set = {p_bg:.5f}")
print(f"\n  {'sample':24s} {'all spec':>10s} " +
      " ".join(f"{'>=' + str(c):>9s}" for c in CUTS) + "   exp n_hp at >=0.35")
for s in samples:
    n = {c: 0 for c in CUTS}; tot = 0
    for c in CH:
        P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
        D = np.load(f"{FEAT}/counts_{s}_chr{c}.npz")
        sil = (P["alt"] == 0) & (P["cov"] >= 15)
        sp = sil & (D["cov"] >= 8) & (D["alt"] >= 2)
        tot += int(sp.sum())
        cov = D["cov"][sp].astype(float); alt = D["alt"][sp].astype(float)
        v = alt / np.maximum(cov, 1)
        for cut in CUTS:
            n[cut] += int((v >= cut).sum())
    exp = n[0.35] * p_bg
    verdict = "POWERLESS" if exp < 10 else ("thin" if exp < 30 else "usable")
    print(f"  {s:24s} {tot:>10,} " + " ".join(f"{n[c]:>9,}" for c in CUTS) +
          f"   {exp:>6.1f}  {verdict}")
print("\n  'exp n_hp' is how many stem>=8 sites you would expect among the CLONAL variants by")
print("  chance. A test whose null expectation is a handful cannot resolve a 1.3x effect.")
