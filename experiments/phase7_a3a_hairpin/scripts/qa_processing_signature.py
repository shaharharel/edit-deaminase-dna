#!/usr/bin/env python
"""Is the PRJNA1006866 control set noisy for a PROCESSING reason? Test it directly.

Section 22 concluded the control-vs-editor split inside one study points at processing
rather than biology. That is a strong claim about provenance and it was inferred from one
number (alt>=1 count). Test it against things processing would move and biology would not:

  STRAND BALANCE    unmarked duplicates or a strand-specific filter skew alt_fwd/alt_rev.
                    Real somatic mutation is strand-balanced at the genome scale.
  SINGLETON SHARE   if the excess sites are alt=1 only, they are sequencing error. If the
                    excess carries alt>=2 as well, error alone does not explain it.
  ALT DEPTH SHAPE   the mean alt count among alt>=1 sites separates "many weak" from
                    "genuinely more mutation".

One chromosome is enough to separate 20x, and chr1 is the largest.
"""
import numpy as np, os
FEAT = os.environ.get("A3A_FEAT", "/mnt/data/a3a/feat")
S = ["A3A-Y130F-clone2", "A3A-Y130F-clone5", "nCas9-clone1", "nCas9-clone2", "Parent",
     "D10A-clone1", "D10A-clone6", "background"]
PRE = {"A3A-Y130F-clone2": "P66-", "A3A-Y130F-clone5": "P66-", "D10A-clone1": "P66-",
       "D10A-clone6": "P66-", "background": "P66-"}
print(f"  {'sample':20s} {'strand f/r':>11} {'alt>=1':>10} {'alt=1 share':>12} "
      f"{'mean alt':>9} {'alt>=2':>10} {'cov':>7}")
for s in S:
    f = f"{FEAT}/counts_{PRE.get(s,'')}{s}_chr1.npz"
    if not os.path.exists(f):
        print(f"  {s:20s} (absent)"); continue
    d = np.load(f)
    cov, alt, af, ar = d["cov"], d["alt"], d["alt_fwd"], d["alt_rev"]
    m = (cov >= 8) & (alt >= 1)
    n1 = int((m & (alt == 1)).sum()); n = int(m.sum()); n2 = int((m & (alt >= 2)).sum())
    sr = af[m].sum() / max(ar[m].sum(), 1)
    print(f"  {s:20s} {sr:11.4f} {n:10,} {100*n1/max(n,1):11.1f}% "
          f"{alt[m].mean():9.3f} {n2:10,} {cov[cov>=8].mean():7.2f}")
print("\n  strand f/r near 1.0 = balanced. A processing artefact skews it; somatic mutation")
print("  does not. alt=1 share near 100% = the excess is sequencing error, not mutation.")
