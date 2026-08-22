#!/usr/bin/env python
"""The specific-site counts in the haA3A output pair up too neatly. Recompute them myself.

    editor Y130G-clone2  90,513   calibrator nCas9-clone1  90,517   -> 4 apart
    editor Y130G-clone1  80,028   calibrator nCas9-clone2  80,053   -> 25 apart
Two independent pairs matching to under 0.03% is the kind of coincidence that is usually a
bug -- the same number printed twice, or crossed labels. Recompute independently from the
counts files instead of trusting the script's own print.

The alt>=1 floors are only ~1% apart for the same pairings (790,444 vs 780,594 and
859,960 vs 854,914), so a genuine burden similarity is plausible -- but "plausible" is not
"checked".
"""
import numpy as np, os
FEAT = os.environ.get("A3A_FEAT", "/data/a3a/feat")
CH = [str(i) for i in range(1, 23)] + ["X"]
PAIRS = [("Y130G-clone2", "nCas9-clone1"), ("Y130G-clone1", "nCas9-clone2"),
         ("Y130G-clone1", "nCas9-clone1")]

def specific(ed, cal, mask="Parent"):
    """editor-specific = jointly eligible, Parent-silent, editor alt>=2, calibrator alt==0."""
    n_ed = n_cal = n_joint = 0
    for c in CH:
        E = np.load(f"{FEAT}/counts_{ed}_chr{c}.npz")
        K = np.load(f"{FEAT}/counts_{cal}_chr{c}.npz")
        P = np.load(f"{FEAT}/counts_{mask}_chr{c}.npz")
        elig = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["alt"] == 0) & (P["cov"] >= 15)
        n_joint += int(elig.sum())
        n_ed  += int((elig & (E["alt"] >= 2) & (K["alt"] == 0)).sum())
        n_cal += int((elig & (K["alt"] >= 2) & (E["alt"] == 0)).sum())
    return n_joint, n_ed, n_cal

print(f"  {'pair':34s} {'jointly elig':>14} {'ed specific':>12} {'cal specific':>13} {'ratio':>7}")
for ed, cal in PAIRS:
    j, a, b = specific(ed, cal)
    print(f"  {ed+' vs '+cal:34s} {j:14,} {a:12,} {b:13,} {a/max(b,1):7.3f}")
print("\n  If the two columns come out DIFFERENT here while the script printed them nearly")
print("  equal, the script has a bug. If they come out nearly equal too, the samples really")
print("  do carry near-identical specific-site burdens and the coincidence is real.")
