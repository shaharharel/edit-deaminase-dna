#!/usr/bin/env python
"""THE THREAT: the editor's mirror read-skew (spread 0.202) is 2.2x the calibrator's (0.090).

Both show the pattern, so it is shared in DIRECTION -- but if it were purely a pipeline
property the magnitudes should match. They do not.

THE HYPOTHESIS THAT WOULD MAKE THIS A REAL THREAT TO THE HEADLINE: hairpin sites sit in
inverted repeats, where alignment is orientation-sensitive, so HAIRPIN SITES may intrinsically
carry more read-orientation skew. The editor set is ENRICHED for hairpins -- that is the
finding -- so its larger skew could be a CONSEQUENCE of that enrichment. But if orientation-
dependent calling is what creates apparent hairpin sites in the first place, the finding is
partly circular.

DECISIVE TEST: compute the skew separately for HAIRPIN and NON-HAIRPIN sites, in BOTH arms.
  - If hairpin sites show the larger skew in BOTH arms -> it is a property of hairpin sites,
    the editor's excess is explained by its composition, and the concern is about whether
    hairpin CALLS are reliable at all.
  - If the editor's skew is larger even WITHIN the non-hairpin stratum -> it is not about
    hairpins, it is about the editor sample.
"""
import numpy as np, sys
FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED, CAL, MASK = sys.argv[1], sys.argv[2], "Parent"
STEM = 6
acc = {(w, h, s): [0, 0] for w in ("ed", "cal") for h in (0, 1) for s in (0, 1)}
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    st, hp = U["strand"], (U["stem"] >= STEM).astype(np.int8)
    elig = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["cov"] >= 15) & (P["alt"] == 0)
    for w, spec, D in (("ed", elig & (E["alt"] >= 2) & (K["alt"] == 0), E),
                       ("cal", elig & (K["alt"] >= 2) & (E["alt"] == 0), K)):
        for h in (0, 1):
            for s in (0, 1):
                m = spec & (hp == h) & (st == s)
                acc[(w, h, s)][0] += int(D["alt_fwd"][m].sum())
                acc[(w, h, s)][1] += int(D["alt_rev"][m].sum())
    del U, E, K, P
print(f"=== mirror read-skew, split by HAIRPIN status (stem>={STEM}) ===")
print(f"  {'arm':5s} {'stratum':13s} {'plus ratio':>11} {'minus ratio':>12} {'spread':>8}")
for w in ("ed", "cal"):
    for h, lab in ((0, "non-hairpin"), (1, "HAIRPIN")):
        r0 = acc[(w,h,0)][0]/max(acc[(w,h,0)][1],1)
        r1 = acc[(w,h,1)][0]/max(acc[(w,h,1)][1],1)
        print(f"  {w:5s} {lab:13s} {r0:11.3f} {r1:12.3f} {abs(r0-r1):8.3f}"
              f"   (n_alt {acc[(w,h,0)][0]+acc[(w,h,0)][1]+acc[(w,h,1)][0]+acc[(w,h,1)][1]:,})")
