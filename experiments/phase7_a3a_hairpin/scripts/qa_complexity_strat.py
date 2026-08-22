#!/usr/bin/env python
"""CHECK 4's UNPERFORMED HALF: does the Y130G signal survive COMPLEXITY stratification?

The standing checklist requires signals to survive GC-decile AND complexity stratification.
GC has been done exhaustively. Complexity never has -- and it is the check most likely to
bite this particular signal, because HAIRPINS ARE INVERTED REPEATS and inverted repeats live
in low-complexity sequence, where mappability and alignment are worst. If the enrichment
concentrates in low-complexity regions it may be an alignment artefact rather than a
deaminase preference.

No raw sequence is stored in the universe files, so k-mer entropy is not available. Instead
use LOCAL TCW SITE DENSITY -- the number of universe sites within +/-500bp of each site.
This is a genuine complexity proxy (TCW sites pack tightly in repeats and spread out in
ordinary sequence), it is computed from `pos` alone, and it is INDEPENDENT of every hairpin
feature, so stratifying on it is not circular the way stratifying on stem length would be.

Same estimator as the GC analysis: per-decile hairpin OR for editor and calibrator, then
Mantel-Haenszel pooling. A shared artefact moves BOTH columns; a real difference does not.
"""
import numpy as np, sys, os
FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED, CAL, MASK = sys.argv[1], sys.argv[2], "Parent"
STEM, WIN = 6, 500

dens_all, hp_all, ed_all, cal_all = [], [], [], []
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    pos, stem = U["pos"].astype(np.int64), U["stem"]
    elig = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["cov"] >= 15) & (P["alt"] == 0)
    # local density: sites within +/-WIN, via searchsorted on the sorted pos array
    lo = np.searchsorted(pos, pos - WIN, side="left")
    hi = np.searchsorted(pos, pos + WIN, side="right")
    dens = (hi - lo).astype(np.int32)
    dens_all.append(dens[elig]); hp_all.append((stem >= STEM)[elig])
    ed_all.append(((E["alt"] >= 2) & (K["alt"] == 0))[elig])
    cal_all.append(((K["alt"] >= 2) & (E["alt"] == 0))[elig])
    del U, E, K, P

d = np.concatenate(dens_all); h = np.concatenate(hp_all)
ed = np.concatenate(ed_all); cal = np.concatenate(cal_all)
print(f"=== {ED} vs {CAL} -- COMPLEXITY (local TCW density in +/-{WIN}bp) stratified ===")
print(f"  jointly eligible {len(d):,}   stem>={STEM} background {h.mean():.5f}")
print(f"  editor specific {int(ed.sum()):,}   calibrator specific {int(cal.sum()):,}")
edges = np.unique(np.percentile(d, np.arange(0, 101, 10)))
b = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, len(edges) - 2)
print(f"\n  {'dec':>3} {'density':>12} {'n_elig':>12} {'p_bg':>8} {'editor':>8} {'n_hp':>7} {'calib':>8} {'n_hp':>7}")
num_e = den_e = num_c = den_c = 0.0
for i in range(len(edges) - 1):
    m = b == i
    if m.sum() < 1000: continue
    pbg = h[m].mean()
    for arr, tag in ((ed, "e"), (cal, "c")):
        s = m & arr; n = int(s.sum()); nh = int((s & h).sum())
        if tag == "e": e_n, e_nh = n, nh
        else: c_n, c_nh = n, nh
    fe = (e_nh / e_n) / pbg if e_n and pbg else float("nan")
    fc = (c_nh / c_n) / pbg if c_n and pbg else float("nan")
    # MH accumulation on the 2x2 (hairpin vs not) x (specific vs background)
    for n_, nh_, which in ((e_n, e_nh, "e"), (c_n, c_nh, "c")):
        if not n_: continue
        a, bb = nh_, n_ - nh_
        cc, dd = pbg * m.sum(), (1 - pbg) * m.sum()
        tot = a + bb + cc + dd
        if which == "e": num_e += a * dd / tot; den_e += bb * cc / tot
        else:            num_c += a * dd / tot; den_c += bb * cc / tot
    print(f"  {i:3d} {edges[i]:6.0f}-{edges[i+1]:<5.0f} {int(m.sum()):12,} {pbg:8.5f} "
          f"{fe:8.3f} {e_nh:7,} {fc:8.3f} {c_nh:7,}")
print(f"\n  MANTEL-HAENSZEL over complexity deciles:  editor {num_e/den_e:.3f}   "
      f"calibrator {num_c/den_c:.3f}   ed-cal {num_e/den_e - num_c/den_c:+.3f}")
print("  A shared complexity artefact moves BOTH columns together. Compare against the")
print("  GC-stratified values (editor 1.357 / calibrator 1.039, ed-cal +0.318) and the")
print("  measured clone-luck floor of +0.075.")
