#!/usr/bin/env python
"""QA of the panel result: replicate on the second clone, and give it a REAL null.

Three things the first run lacked:
  1. NO RANDOM BASELINE AT 0.01% -- the depth where the effect is largest and n smallest.
     That is the non-negotiable check, missing exactly where it matters most.
  2. NO REPLICATION. Y130G-clone1 is an independent clone and the whole Y130G story rests
     on the two clones agreeing.
  3. NO PERMUTATION NULL. A single shuffled draw is not a null distribution; with n=12 in a
     cell, one draw says nothing. Use 200 permutations and report the null's spread.
"""
import numpy as np, os, sys, time
from sklearn.ensemble import HistGradientBoostingClassifier
FEAT="/data/a3a/feat"; CH=[str(i) for i in range(1,23)]+["X"]
ED=sys.argv[1]; CAL=sys.argv[2]; MASK="Parent"; SEED=0; NPERM=200
FEATS=("stem","loop","kpos","gc_pairs","hp_score","local_gc")
t0=time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}",flush=True)
d=np.load(f"{FEAT}/a3a_trainset_v5.npz",allow_pickle=True)
Xtr=np.column_stack([d[f] for f in FEATS]).astype(np.float32); ytr=d["y"].astype(np.int8)
m=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.1,early_stopping=True,
                                 validation_fraction=0.1,random_state=SEED).fit(Xtr,ytr)
log("PCAWG model fitted")
sc,ed_s,cal_s=[],[],[]
for c in CH:
    U=np.load(f"{FEAT}/universe_chr{c}.npz"); E=np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K=np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz"); P=np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    el=(E["cov"]>=8)&(K["cov"]>=8)&(P["cov"]>=15)&(P["alt"]==0)
    X=np.column_stack([U[f][el] for f in FEATS]).astype(np.float32)
    sc.append(m.predict_proba(X)[:,1].astype(np.float32))
    ed_s.append(((E["alt"]>=2)&(K["alt"]==0))[el]); cal_s.append(((K["alt"]>=2)&(E["alt"]==0))[el])
    del U,E,K,P,X
sc=np.concatenate(sc); ed=np.concatenate(ed_s); cal=np.concatenate(cal_s); N=len(sc)
log(f"{ED}: {N:,} sites, ed-specific {int(ed.sum()):,}, cal-specific {int(cal.sum()):,}")
order=np.argsort(-sc); rng=np.random.default_rng(SEED)
edi=np.flatnonzero(ed); cali=np.flatnonzero(cal)
print(f"\n  {'panel':>7} {'n_ed':>6} {'ed x':>7} {'n_cal':>6} {'cal x':>7} {'gap':>8} "
      f"{'null gap mean':>14} {'null sd':>8} {'z':>6}")
for k in (0.01,0.10,1.00,5.00):
    n=max(1,int(N*k/100)); top=order[:n]
    tm=np.zeros(N,bool); tm[top]=True
    ne=int(tm[edi].sum()); nc=int(tm[cali].sum())
    ex=(ne/len(edi))/(k/100); cx=(nc/len(cali))/(k/100); gap=ex-cx
    gaps=np.empty(NPERM)
    for i in range(NPERM):
        r=rng.choice(N,size=n,replace=False)
        rm=np.zeros(N,bool); rm[r]=True
        gaps[i]=((rm[edi].sum()/len(edi))-(rm[cali].sum()/len(cali)))/(k/100)
    z=(gap-gaps.mean())/max(gaps.std(ddof=1),1e-9)
    print(f"  {k:6.2f}% {ne:6d} {ex:7.3f} {nc:6d} {cx:7.3f} {gap:+8.3f} "
          f"{gaps.mean():+14.3f} {gaps.std(ddof=1):8.3f} {z:6.1f}")
print(f"\n  null = {NPERM} permutations of the RANKING at each panel size, gap recomputed each time.")
