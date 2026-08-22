#!/usr/bin/env python
"""IS THE EDITOR DATA LEARNABLE AT ALL? Train and test ON THE EDITOR ITSELF.

Everything so far trained on PCAWG and asked whether it transfers. That question is now
answered and the answer is that the cancer training added nothing over raw `stem`. The
DIFFERENT question -- never asked -- is whether the editor's own sites are learnable from
sequence/structure at all, when the model is allowed to train on them.

DESIGN
  positives  editor-specific sites (editor alt>=2, calibrator alt==0, Parent-silent)
  negatives  sampled from eligible sites the editor did NOT call, MATCHED ON TRINUCLEOTIDE
             AND STRAND, which is non-negotiable in this project
  split      held-out CHROMOSOME, never random
  model      HistGradientBoosting on the 6 hairpin/structure features
  endpoint   TAIL ENRICHMENT with a RANDOM baseline beside it, not AUROC

AND THE CONTROL THAT DECIDES WHETHER ANY OF IT MEANS "EDITOR": run the identical pipeline on
the CALIBRATOR's own specific sites. nCas9 is deaminase-free, so whatever the model learns
there is endogenous. If editor and calibrator are equally learnable, the model is learning
somatic mutation structure, NOT editor activity.
"""
import numpy as np, sys, time
from sklearn.ensemble import HistGradientBoostingClassifier
FEAT="/data/a3a/feat"; CH=[str(i) for i in range(1,23)]+["X"]
ED=sys.argv[1]; CAL=sys.argv[2]; MASK="Parent"; SEED=0; NEG_RATIO=10
FEATS=("stem","loop","kpos","gc_pairs","hp_score","local_gc")
t0=time.time(); rng=np.random.default_rng(SEED)
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}",flush=True)

def build(which):
    """which='editor' -> editor-specific positives; 'calibrator' -> calibrator-specific."""
    Xs,ys,chs=[],[],[]
    for c in CH:
        U=np.load(f"{FEAT}/universe_chr{c}.npz"); E=np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
        K=np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz"); P=np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
        el=(E["cov"]>=8)&(K["cov"]>=8)&(P["cov"]>=15)&(P["alt"]==0)
        pos = ((E["alt"]>=2)&(K["alt"]==0))&el if which=="editor" else ((K["alt"]>=2)&(E["alt"]==0))&el
        neg_pool = el & ~pos & (E["alt"]==0) & (K["alt"]==0)
        tri,st=U["tri"],U["strand"]
        keep_p=np.flatnonzero(pos); sel_n=[]
        for t in (0,1):
            for s in (0,1):
                need=int(((tri[keep_p]==t)&(st[keep_p]==s)).sum())*NEG_RATIO
                cand=np.flatnonzero(neg_pool&(tri==t)&(st==s))
                if need and len(cand): sel_n.append(rng.choice(cand,size=min(need,len(cand)),replace=False))
        keep_n=np.concatenate(sel_n) if sel_n else np.array([],dtype=int)
        idx=np.concatenate([keep_p,keep_n])
        Xs.append(np.column_stack([U[f][idx] for f in FEATS]).astype(np.float32))
        ys.append(np.concatenate([np.ones(len(keep_p),np.int8),np.zeros(len(keep_n),np.int8)]))
        chs.append(np.full(len(idx),c,dtype=object))
        del U,E,K,P
    return np.vstack(Xs),np.concatenate(ys),np.concatenate(chs)

def run(which):
    X,y,ch=build(which)
    log(f"{which}: n={len(y):,} pos={int(y.sum()):,} base_rate={y.mean():.4f} ceiling={1/y.mean():.2f}x")
    # leak check: trinuc+strand must match between classes
    log(f"  {which} leak check -- base rate {y.mean():.4f}; negatives matched on trinuc AND strand by construction")
    chroms=sorted(set(ch.tolist()),key=lambda c:(len(c),c)); folds=[chroms[i::5] for i in range(5)]
    oof=np.zeros(len(y),np.float32)
    for hold in folds:
        te=np.isin(ch,hold); tr=~te
        if te.sum()==0 or y[tr].sum()==0: continue
        mm=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.1,early_stopping=True,
                                          validation_fraction=0.1,random_state=SEED).fit(X[tr],y[tr])
        oof[te]=mm.predict_proba(X[te])[:,1]
    base=y.mean()
    print(f"  {which:11s} {'topK':>7} {'enr':>8} {'RANDOM':>8} {'n':>8}")
    for k in (0.1,1.0,5.0):
        n=max(1,int(len(y)*k/100))
        top=np.argsort(-oof)[:n]; e=y[top].mean()/base
        r=y[rng.permutation(len(y))[:n]].mean()/base
        print(f"  {'':11s} {k:6.1f}% {e:8.3f} {r:8.3f} {n:8,}")
    return

for w in ("editor","calibrator"):
    run(w)
print("\n  If EDITOR and CALIBRATOR are equally learnable, the model is learning somatic")
print("  mutation structure, not editor activity. The calibrator is deaminase-free.")
