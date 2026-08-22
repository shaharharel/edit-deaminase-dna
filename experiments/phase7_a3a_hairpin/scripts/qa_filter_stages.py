#!/usr/bin/env python
"""WHERE DOES THE SPECIFICITY FILTER LOSE EACH ARM'S SITES? (§27's blocking question)

At alt>=2 the A3A-Y130F clones hold 72% of what nCas9 holds. After the specificity filter
they hold 22%. The gap TRIPLES in the step that defines the analysis set, and the obvious
mechanism -- the calibrator veto -- cannot account for it (the calibrator's own alt>=2 rate
is 0.12% of eligible sites).

Decompose the filter stage by stage, separately per arm, and find the step that does it:
    S0  editor alt>=2, editor cov>=8                     (raw candidate set)
    S1  + calibrator cov>=8                              (joint eligibility)
    S2  + Parent cov>=15                                 (mask has power here)
    S3  + Parent alt==0                                  (germline mask applied)
    S4  + calibrator alt==0                              (the calibrator veto)  = SPECIFIC
Report survival at every step for both arms. Whichever step diverges IS the explanation.
"""
import numpy as np, os
FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ARMS = [("P66-A3A-Y130F-clone2", "nCas9-clone1"),
        ("P66-A3A-Y130F-clone5", "nCas9-clone1"),
        ("nCas9-clone2",         "nCas9-clone1"),
        ("P66-A3A-Y130F-clone2", "P66-D10A-clone6")]

for ed, cal in ARMS:
    s = np.zeros(5, dtype=np.int64)
    for c in CH:
        E = np.load(f"{FEAT}/counts_{ed}_chr{c}.npz")
        K = np.load(f"{FEAT}/counts_{cal}_chr{c}.npz")
        P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
        m = (E["cov"] >= 8) & (E["alt"] >= 2);            s[0] += int(m.sum())
        m &= (K["cov"] >= 8);                             s[1] += int(m.sum())
        m &= (P["cov"] >= 15);                            s[2] += int(m.sum())
        m &= (P["alt"] == 0);                             s[3] += int(m.sum())
        m &= (K["alt"] == 0);                             s[4] += int(m.sum())
    lab = ["S0 ed alt>=2", "S1 +cal cov>=8", "S2 +Par cov>=15", "S3 +Par alt==0", "S4 +cal veto"]
    print(f"\n=== {ed}  vs  {cal} ===")
    print(f"  {'stage':18s} {'surviving':>12} {'% of S0':>9} {'lost at this step':>18}")
    for i in range(5):
        drop = (s[i-1] - s[i]) if i else 0
        print(f"  {lab[i]:18s} {s[i]:12,} {100*s[i]/max(s[0],1):8.1f}% "
              f"{(f'{drop:,}' if i else '--'):>18}")
