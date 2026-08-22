#!/usr/bin/env python
"""s14 v2 -- GENOMIC-SCALE CANCER SEARCH, OUT-OF-FOLD. v1 was 90.5% in-sample.

v1 fitted ONE model on all of v5 and scored the whole genome, then joined mutations from the
SAME donors. 83,977 of the 92,823 scored positives were training positives. Its 7.337x was
memorisation, and comparing it to the editor's 1.748x was unfair -- the editor arm was scored
by a model trained on PCAWG that had never seen a HEK293T site, i.e. genuinely out-of-sample.

FIX AT THE GENERATOR: held-out CHROMOSOME folds, the same rule every other model in this
project obeys. Train on the trainset rows from 4/5 of chromosomes, score the universe sites
of the held-out fifth, and assemble genome-wide scores that are OUT-OF-FOLD everywhere.
Then the cancer number and the editor number are both out-of-sample and comparable.
"""
import numpy as np, time, sys
from sklearn.ensemble import HistGradientBoostingClassifier
FEAT="/data/a3a/feat"; PC="/data/a3a/pcawg"
CH=[str(i) for i in range(1,23)]+["X"]
FEATS=("stem","loop","kpos","gc_pairs","hp_score","local_gc"); SEED=0
t0=time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}",flush=True)

d=np.load(f"{FEAT}/a3a_trainset_v5.npz",allow_pickle=True)
Xtr=np.column_stack([d[f] for f in FEATS]).astype(np.float32)
ytr=d["y"].astype(np.int8); ctr=d["chrom"]
train_donors=set(np.unique(d["donor"]).tolist())
chroms=sorted(set(ctr.tolist()),key=lambda c:(len(c),c))
folds=[chroms[i::5] for i in range(5)]
log(f"5 folds by chromosome: {[len(f) for f in folds]} chroms each")

s=np.load(f"{PC}/snvs.npz",allow_pickle=True)
keep=np.isin(s["donor"],list(train_donors))
mu_chrom=s["chrom"][keep]; mu_pos=s["pos"][keep].astype(np.int64)

scores=[]; ispos=[]; tot_mut=0
for fi,hold in enumerate(folds):
    tr=~np.isin(ctr,hold)
    m=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.1,early_stopping=True,
                                     validation_fraction=0.1,random_state=SEED).fit(Xtr[tr],ytr[tr])
    for c in hold:
        U=np.load(f"{FEAT}/universe_chr{c}.npz")
        upos=U["pos1"] if "pos1" in U.files else U["pos"].astype(np.int64)+1
        X=np.column_stack([U[f] for f in FEATS]).astype(np.float32)
        scores.append(m.predict_proba(X)[:,1].astype(np.float32))
        mp=np.sort(mu_pos[mu_chrom==c])
        i=np.searchsorted(upos,mp); i=np.clip(i,0,len(upos)-1)
        hit=np.zeros(len(upos),bool); good=upos[i]==mp
        hit[i[good]]=True; tot_mut+=int(hit.sum())
        ispos.append(hit); del U,X
    log(f"fold {fi}: held out {len(hold)} chroms, model never saw them")
sc=np.concatenate(scores); pos=np.concatenate(ispos); N=len(sc)
br=pos.mean()
log(f"universe {N:,} sites, OUT-OF-FOLD scored; mutations on universe {tot_mut:,}")
print(f"\n  GENOMIC base rate {br:.6f} = 1 in {1/br:,.0f}")
print(f"  EVERY score is out-of-fold: the model that scored a site never trained on its chromosome.")
rng=np.random.default_rng(SEED); order=np.argsort(-sc)
print(f"\n  {'panel':>7} {'sites':>13} {'found':>10} {'capture':>9} {'enrich':>8} {'RANDOM':>8}")
for k in (0.01,0.1,1.0,5.0,10.0):
    n=max(1,int(N*k/100)); c_=int(pos[order[:n]].sum())
    r_=int(pos[rng.permutation(N)[:n]].sum())
    print(f"  {k:6.2f}% {n:13,} {c_:10,} {100*c_/tot_mut:8.2f}% "
          f"{(100*c_/tot_mut)/k:8.3f} {(100*r_/tot_mut)/k:8.3f}")
print("\n  compare to the EDITOR at genomic scale (also out-of-sample):")
print("    top 0.1%  1.748x    top 1.0%  1.146x    top 5.0%  0.985x")
