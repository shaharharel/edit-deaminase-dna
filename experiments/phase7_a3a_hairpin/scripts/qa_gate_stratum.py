#!/usr/bin/env python
"""GATE A0 JUDGES alt>=1. THE ANALYSIS USES alt>=2. Those may not be the same sample.

chr1 showed the disqualified samples' excess is overwhelmingly SINGLETONS:
    alt=1 share   clean 70-80%   D10A-clone1 94.6%   D10A-clone6 94.0%   background 88.5%
    mean alt      clean 3.6-4.0  D10A 1.49 / 1.67    background 1.34
and their alt>=2 counts are NOT similarly inflated:
    alt>=2   clean 15,430-27,516   D10A 33,321 / 31,137   background 176,912
D10A looks proportional to its greater depth; background does not.

If a sample's defect lives entirely in the alt=1 stratum, and s7c_editor's analysis set is
alt>=2, then the gate is rejecting on a stratum the analysis never reads. That is this
project's bug family pointed at my own gate: a threshold correct in one place applied where
a different one governs.

Recompute the depth-adjusted residual on alt>=2 -- the stratum the analysis actually uses --
over ALL 23 chromosomes, and see which disqualifications survive.
"""
import numpy as np, os
FEAT = os.environ.get("A3A_FEAT", "/mnt/data/a3a/feat")
CH = [str(i) for i in range(1, 23)] + ["X"]
SAMP = ["P66-A3A-Y130F-clone2", "P66-A3A-Y130F-clone5", "nCas9-clone1", "nCas9-clone2",
        "Parent", "P66-D10A-clone1", "P66-D10A-clone6", "P66-background"]
ROLE = {"nCas9-clone1": "control", "nCas9-clone2": "control", "Parent": "control",
        "P66-D10A-clone1": "control", "P66-D10A-clone6": "control",
        "P66-background": "control", "P66-A3A-Y130F-clone2": "EDITOR",
        "P66-A3A-Y130F-clone5": "EDITOR"}

res = {}
for s in SAMP:
    n1 = n2 = 0; covs = 0.0; nsite = 0
    for c in CH:
        f = f"{FEAT}/counts_{s}_chr{c}.npz"
        if not os.path.exists(f): n1 = -1; break
        d = np.load(f); cov, alt = d["cov"], d["alt"]
        m = cov >= 8
        n1 += int((m & (alt >= 1)).sum()); n2 += int((m & (alt >= 2)).sum())
        covs += float(cov[m].sum()); nsite += int(m.sum())
    if n1 < 0: print(f"  {s}: incomplete"); continue
    res[s] = (n1, n2, covs / max(nsite, 1))

# fit the depth trend on the SAME five samples the v3 gate used as sound
SOUND = ["nCas9-clone1", "nCas9-clone2", "Parent",
         "P66-A3A-Y130F-clone2", "P66-A3A-Y130F-clone5"]
for lab, idx in (("alt>=1  (what the gate judges)", 0), ("alt>=2  (what the ANALYSIS uses)", 1)):
    c = np.array([res[s][2] for s in SOUND]); n = np.array([res[s][idx] for s in SOUND], float)
    b, a = np.polyfit(c, np.log10(n), 1)
    r = np.corrcoef(c, np.log10(n))[0, 1]
    print(f"\n=== {lab}   trend fitted on the 5 sound samples, r={r:.4f} ===")
    print(f"  {'sample':22s} {'role':8s} {'cov':>6} {'observed':>12} {'expected':>12} {'resid':>8}")
    for s in SAMP:
        if s not in res: continue
        n_o, cv = res[s][idx], res[s][2]
        exp = 10 ** (a + b * cv)
        flag = "  <-- FAILS 2.0x" if n_o / exp > 2.0 and ROLE[s] == "control" else ""
        print(f"  {s:22s} {ROLE[s]:8s} {cv:6.2f} {n_o:12,} {exp:12,.0f} {n_o/exp:7.2f}x{flag}")
