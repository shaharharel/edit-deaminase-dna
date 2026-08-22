#!/usr/bin/env python
"""CHECK 1 ON THE LEARNABILITY TEST -- which I asserted instead of measuring.

s13_learnable.py printed "negatives matched on trinuc AND strand by construction". I did not
VERIFY it. "By construction" is exactly what I said about the atomic write (which silently
mis-named 23 files) and about the peer list (which silently fell back on node B). It is not
evidence.

There is a specific way it could be false. The matcher does
    need = (#positives in this tri x strand stratum) * NEG_RATIO
    take = min(need, len(candidates))
so IF ANY STRATUM RAN SHORT OF CANDIDATES, the ratio silently drops for that stratum only and
the trinucleotide/strand balance between classes breaks -- with no error and no warning.

Measure it: rebuild the same sets and compare the tri and strand distributions class by class,
plus the realised negative:positive ratio per stratum.
"""
import numpy as np, sys
FEAT="/data/a3a/feat"; CH=[str(i) for i in range(1,23)]+["X"]
ED=sys.argv[1]; CAL=sys.argv[2]; MASK="Parent"; SEED=0; NEG_RATIO=10
rng=np.random.default_rng(SEED)

for which in ("editor","calibrator"):
    npos=np.zeros((2,2),dtype=np.int64); nneg=np.zeros((2,2),dtype=np.int64)
    short=[]
    for c in CH:
        U=np.load(f"{FEAT}/universe_chr{c}.npz"); E=np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
        K=np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz"); P=np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
        el=(E["cov"]>=8)&(K["cov"]>=8)&(P["cov"]>=15)&(P["alt"]==0)
        pos=(((E["alt"]>=2)&(K["alt"]==0)) if which=="editor" else ((K["alt"]>=2)&(E["alt"]==0)))&el
        neg_pool=el&~pos&(E["alt"]==0)&(K["alt"]==0)
        tri,st=U["tri"],U["strand"]
        kp=np.flatnonzero(pos)
        for t in (0,1):
            for s in (0,1):
                n_p=int(((tri[kp]==t)&(st[kp]==s)).sum()); npos[t,s]+=n_p
                need=n_p*NEG_RATIO
                cand=np.flatnonzero(neg_pool&(tri==t)&(st==s))
                take=min(need,len(cand)); nneg[t,s]+=take
                if need and take<need: short.append((c,t,s,need,take))
        del U,E,K,P
    print(f"\n=== {which} ===")
    tp,tn=npos.sum(),nneg.sum()
    print(f"  {'tri':>4} {'strand':>7} {'positives':>11} {'negatives':>11} {'realised ratio':>15} "
          f"{'pos share':>10} {'neg share':>10}")
    for t,tl in ((0,"TCA"),(1,"TCT")):
        for s in (0,1):
            p_,n_=npos[t,s],nneg[t,s]
            print(f"  {tl:>4} {s:>7} {p_:11,} {n_:11,} {n_/max(p_,1):15.3f} "
                  f"{p_/tp:10.4f} {n_/tn:10.4f}")
    d=max(abs(npos[t,s]/tp - nneg[t,s]/tn) for t in (0,1) for s in (0,1))
    print(f"  MAX |pos share - neg share| across the four strata = {d:.5f} "
          f"{'OK' if d<0.005 else '*** MATCHING IS BROKEN ***'}")
    print(f"  strata that ran SHORT of candidates: {len(short)}"
          f"{' -- ' + str(short[:3]) if short else ' (none, so the 10:1 ratio held everywhere)'}")
