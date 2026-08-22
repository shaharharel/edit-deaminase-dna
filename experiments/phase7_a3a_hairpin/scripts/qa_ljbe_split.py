#!/usr/bin/env python
"""clone3 has 1,299,830 specific sites and clone5 has 201,637 -- a 6.4x difference between two
clones of the SAME arm. qa_2way_ljbe.py reported the private fraction for clone3 only, and
clone3 supplies 1,281,842 of the 1,299,830 private sites, so "private = 1.369x" is essentially
clone3's number. If clone3 simply has a higher noise floor, that is what is being measured.

Splits the private fraction by clone, and checks each clone's raw noise floor with the same
criterion that disqualified D10A-clone1.
"""
import numpy as np

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
A, B = "P66-Lj-BE-clone3", "P66-Lj-BE-clone5"
STEMS = [6, 7, 8]

pa = {s: [0, 0] for s in STEMS}
pb = {s: [0, 0] for s in STEMS}
bg = {s: [0, 0] for s in STEMS}
nf = {A: [0, 0], B: [0, 0]}

for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    Da = np.load(f"{FEAT}/counts_{A}_chr{c}.npz")
    Db = np.load(f"{FEAT}/counts_{B}_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15)
    el = sil & (Da["cov"] >= 8) & (Db["cov"] >= 8)
    sa = el & (Da["alt"] >= 2); sb = el & (Db["alt"] >= 2)
    onlya = sa & ~sb; onlyb = sb & ~sa
    for s in STEMS:
        h = stem >= s
        bg[s][0] += int((el & h).sum()); bg[s][1] += int(el.sum())
        pa[s][0] += int((onlya & h).sum()); pa[s][1] += int(onlya.sum())
        pb[s][0] += int((onlyb & h).sum()); pb[s][1] += int(onlyb.sum())
    for nm, D in ((A, Da), (B, Db)):
        m = (D["alt"] >= 1) & (D["cov"] >= 8)
        nf[nm][0] += int(m.sum())
        v = D["alt"][m] / np.maximum(D["cov"][m], 1)
        nf[nm][1] += int((v < 0.05).sum())

print("=== raw noise floor, the criterion that disqualified D10A-clone1 ===")
for nm in (A, B):
    n, lo = nf[nm]
    print(f"  {nm:22s} alt>=1 {n:>11,}   {100*lo/max(n,1):5.1f}% below VAF 0.05")
print("  healthy band measured on node A: 780k-1.14M sites, 51.8-67.2% sub-0.05")

print("\n=== private fraction SPLIT BY CLONE (qa_2way_ljbe reported clone3 only) ===")
print(f"  {'stem':>5} {'p_bg':>9} {'clone3-only':>12} {'n_hp':>8} {'clone5-only':>12} {'n_hp':>7}")
for s in STEMS:
    p = bg[s][0] / max(bg[s][1], 1)
    ea = (pa[s][0] / pa[s][1]) / p if pa[s][1] else float("nan")
    eb = (pb[s][0] / pb[s][1]) / p if pb[s][1] else float("nan")
    print(f"  {s:>5} {p:>9.5f} {ea:>11.3f}x {pa[s][0]:>8,} {eb:>11.3f}x {pb[s][0]:>7,}")
print(f"  n: clone3-only {pa[8][1]:,}   clone5-only {pb[8][1]:,}")
print("  eA3A private for comparison: 1.371x at stem 8. Calibrator range 1.28-1.63x.")
