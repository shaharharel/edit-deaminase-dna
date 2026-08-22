#!/usr/bin/env python
"""Why does Lj-BE-clone3 call 6.4x more specific sites than clone5?

The comparison is already apples-to-apples: both counts come from the SAME jointly-eligible
set (Parent-silent, cov>=8 in BOTH clones), so it is not an eligibility difference. And the
direction is wrong for a coverage explanation -- clone5 is DEEPER (41.8x vs 38.0x), so at a
fixed alt>=2 rule it should call MORE sites, not 6.4x fewer.

Their raw noise floors differ 2.7x, which leaves a factor of ~2.4 unexplained.

The discriminator is the VAF spectrum of the specific sites. If clone3's excess sits at low
VAF it is background; if it is spread across the spectrum it is a real difference in mutation
load between two clones of the same arm.
"""
import numpy as np

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
A, B = "P66-Lj-BE-clone3", "P66-Lj-BE-clone5"
EDGES = [0.0, 0.05, 0.10, 0.15, 0.25, 0.35, 1.01]

prof = {A: np.zeros(len(EDGES) - 1, np.int64), B: np.zeros(len(EDGES) - 1, np.int64)}
tot = {A: 0, B: 0}
elig = 0
for c in CH:
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    Da = np.load(f"{FEAT}/counts_{A}_chr{c}.npz")
    Db = np.load(f"{FEAT}/counts_{B}_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15)
    el = sil & (Da["cov"] >= 8) & (Db["cov"] >= 8)     # identical set for both
    elig += int(el.sum())
    for nm, D in ((A, Da), (B, Db)):
        sp = el & (D["alt"] >= 2)
        tot[nm] += int(sp.sum())
        cov = D["cov"][sp].astype(float); alt = D["alt"][sp].astype(float)
        v = alt / np.maximum(cov, 1)
        for i in range(len(EDGES) - 1):
            prof[nm][i] += int(((v >= EDGES[i]) & (v < EDGES[i + 1])).sum())

print(f"jointly eligible sites (identical for both): {elig:,}")
hdr = "  ".join(f"{EDGES[i]:.2f}-{EDGES[i+1]:.2f}" for i in range(len(EDGES) - 1))
print(f"\n  {'clone':22s} {'n_spec':>10s}   {hdr}")
for nm in (A, B):
    p = prof[nm]; t = max(p.sum(), 1)
    print(f"  {nm:22s} {tot[nm]:>10,}   " + "  ".join(f"{100*x/t:8.1f}%" for x in p))
print(f"\n  ABSOLUTE counts per VAF bin (the ratio is what matters):")
print(f"  {'bin':>12s} {'clone3':>12s} {'clone5':>12s} {'ratio':>8s}")
for i in range(len(EDGES) - 1):
    a, b = int(prof[A][i]), int(prof[B][i])
    r = a / b if b else float("nan")
    print(f"  {EDGES[i]:.2f}-{EDGES[i+1]:.2f} {a:>12,} {b:>12,} {r:>7.2f}x")
print("\n  If the excess is confined to the lowest bins it is background. If the ratio is flat")
print("  across the spectrum, the two clones genuinely differ in mutation load.")
