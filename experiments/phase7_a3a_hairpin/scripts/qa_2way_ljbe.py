#!/usr/bin/env python
"""TWO-clone recurrence test on the Lj-BE arm, the moment clone5's counts land.

This is the instrument that downgraded the A3A-Y130F arm. There, two independent clones shared
3,000 editor-specific sites against 1.8 expected by chance (1642x), and splitting the arm on
that sharing showed the marginal stem-8 signal came ENTIRELY from the shared fraction
(1.326x shared vs 1.026x private). One clone could never have revealed it.

The Lj-BE single-clone preview said lamprey CDA1 matches eA3A, which under the pre-registration
is the "substrate, not enzyme" branch. That preview is exactly the kind of claim this test
exists to check, and it needs only TWO clones -- so there is no reason to wait for clone12,
which is still ~3 h out.

PRE-REGISTERED, before running:
  - if the Lj-BE signal sits in the PRIVATE fraction, it behaves like eA3A (private 1.371x vs
    all-three 1.110x at stem 8) and the preview stands;
  - if it sits in the SHARED fraction, it is the A3A-Y130F artefact again and the preview was
    misleading.
Either way the numbers are quoted against the deaminase-free calibrator's 1.28-1.63x, never
against 1.0, and n is reported with every figure.
"""
import numpy as np

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
A, B = "P66-Lj-BE-clone3", "P66-Lj-BE-clone5"
STEMS = [6, 7, 8]

ov = {"a": 0, "b": 0, "both": 0, "elig": 0}
hp = {s: {"shared": [0, 0], "private_a": [0, 0]} for s in STEMS}
bg = {s: [0, 0] for s in STEMS}

for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    Da = np.load(f"{FEAT}/counts_{A}_chr{c}.npz")
    Db = np.load(f"{FEAT}/counts_{B}_chr{c}.npz")
    silent = (P["alt"] == 0) & (P["cov"] >= 15)
    el = silent & (Da["cov"] >= 8) & (Db["cov"] >= 8)      # JOINTLY eligible
    sa = el & (Da["alt"] >= 2)
    sb = el & (Db["alt"] >= 2)
    both = sa & sb
    ov["elig"] += int(el.sum()); ov["a"] += int(sa.sum())
    ov["b"] += int(sb.sum());   ov["both"] += int(both.sum())
    for s in STEMS:
        h = stem >= s
        bg[s][0] += int((el & h).sum()); bg[s][1] += int(el.sum())
        hp[s]["shared"][0] += int((both & h).sum());  hp[s]["shared"][1] += int(both.sum())
        pv = sa & ~sb
        hp[s]["private_a"][0] += int((pv & h).sum()); hp[s]["private_a"][1] += int(pv.sum())

exp = ov["a"] * ov["b"] / max(ov["elig"], 1)
print(f"=== cross-clone recurrence, {A} vs {B} ===")
print(f"  jointly eligible   {ov['elig']:>12,}")
print(f"  {A:18s} {ov['a']:>12,}")
print(f"  {B:18s} {ov['b']:>12,}")
print(f"  shared by BOTH     {ov['both']:>12,}")
print(f"  expected by chance {exp:>12,.1f}")
print(f"  OBS/EXP            {ov['both']/max(exp,1e-9):>12.1f}x    "
      f"(A3A-Y130F was 1642x on 3,000 shared sites)")

print(f"\n=== hairpin enrichment, SHARED vs PRIVATE (the split that downgraded A3A-Y130F) ===")
print(f"  {'stem':>5} {'p_bg':>9} {'shared':>9} {'n_hp':>7} {'private':>9} {'n_hp':>7}")
for s in STEMS:
    pbg = bg[s][0] / max(bg[s][1], 1)
    sh_h, sh_n = hp[s]["shared"]
    pv_h, pv_n = hp[s]["private_a"]
    e_sh = (sh_h / sh_n) / pbg if sh_n else float("nan")
    e_pv = (pv_h / pv_n) / pbg if pv_n else float("nan")
    print(f"  {s:>5} {pbg:>9.5f} {e_sh:>8.3f}x {sh_h:>7,} {e_pv:>8.3f}x {pv_h:>7,}")
print(f"\n  shared n={ov['both']:,}  private(clone3) n={ov['a']-ov['both']:,}")
print(f"  Read against the deaminase-free calibrator 1.28-1.63x, never against 1.0.")
print(f"  eA3A for comparison: private 1.371x vs all-three 1.110x at stem 8 -- signal in the")
print(f"  PRIVATE fraction. A3A-Y130F was the opposite and was downgraded for it.")
