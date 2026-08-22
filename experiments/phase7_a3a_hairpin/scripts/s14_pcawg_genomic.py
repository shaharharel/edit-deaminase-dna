#!/usr/bin/env python
"""GENOMIC-SCALE SEARCH IN CANCER -- the gap I just admitted to the user and never measured.

I have been quoting 5.280x for the cancer model. That is enrichment over a base rate of
0.0909 THAT I CONSTRUCTED by sampling 10 matched negatives per positive. It measures
DISCRIMINATION: editor... sorry, MUTATED site vs MATCHED non-mutated site.

It does NOT measure what a panel would need: score every TCW site in the genome, take the
top K%, and ask what fraction of real PCAWG mutations land there. Today's HEK293T work showed
those two scales differ by 234x, so the question is not academic -- and I told the user this
number was unmeasured. Measuring it.

DESIGN
  score      every site in the hg19 TCW universe with the SAME PCAWG-trained model
  positives  real PCAWG APOBEC mutations from snvs.npz, restricted to the SAME donor set the
             model was trained on, mapped onto the universe by (chrom, pos)
  endpoint   capture at top K%, with the RANDOM baseline beside it
  and the ceiling recomputed for the GENOMIC base rate, which is not 11.00x

COORDINATE CONVENTION -- bug 3 recurred 22 times in this project, so it is handled explicitly:
the universe stores `pos1` (1-based) alongside 0-based `pos`; snvs.npz `pos` is 1-based.
Join on pos1 and VERIFY the join rate, refusing to report if it is implausible.
"""
import numpy as np, time, sys
from sklearn.ensemble import HistGradientBoostingClassifier
FEAT="/data/a3a/feat"; PC="/data/a3a/pcawg"
CH=[str(i) for i in range(1,23)]+["X"]
FEATS=("stem","loop","kpos","gc_pairs","hp_score","local_gc"); SEED=0
t0=time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}",flush=True)

d=np.load(f"{FEAT}/a3a_trainset_v5.npz",allow_pickle=True)
Xtr=np.column_stack([d[f] for f in FEATS]).astype(np.float32); ytr=d["y"].astype(np.int8)
train_donors=set(np.unique(d["donor"]).tolist())
m=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.1,early_stopping=True,
                                 validation_fraction=0.1,random_state=SEED).fit(Xtr,ytr)
log(f"model fitted on v5; {len(train_donors)} training donors")

s=np.load(f"{PC}/snvs.npz",allow_pickle=True)
keep=np.isin(s["donor"],list(train_donors))
mut_chrom=s["chrom"][keep]; mut_pos=s["pos"][keep].astype(np.int64)
log(f"PCAWG mutations from those donors: {keep.sum():,} of {len(keep):,}")

tot_sites=0; tot_mut=0; scores=[]; ispos=[]
for c in CH:
    U=np.load(f"{FEAT}/universe_chr{c}.npz")
    upos=U["pos1"] if "pos1" in U.files else U["pos"].astype(np.int64)+1
    X=np.column_stack([U[f] for f in FEATS]).astype(np.float32)
    sc=m.predict_proba(X)[:,1].astype(np.float32)
    mp=np.sort(mut_pos[mut_chrom==c])
    idx=np.searchsorted(upos,mp)
    ok=(idx<len(upos)); idx=idx[ok]
    hit=np.zeros(len(upos),bool)
    good=upos[idx]==mp[ok]
    hit[idx[good]]=True
    tot_sites+=len(upos); tot_mut+=int(hit.sum())
    scores.append(sc); ispos.append(hit)
    del U,X
sc=np.concatenate(scores); pos=np.concatenate(ispos)
join=tot_mut/max(keep.sum(),1)
log(f"universe {tot_sites:,} sites; PCAWG mutations landing on it {tot_mut:,} "
    f"({100*join:.1f}% of the donor set)")
if tot_mut < 1000:
    sys.exit("*** JOIN FAILED -- refusing to report, check the coordinate convention ***")
br=pos.mean()
print(f"\n  GENOMIC base rate {br:.6f} = 1 in {1/br:,.0f}   ceiling {1/br:,.1f}x")
print(f"  (the trainset base rate was 0.0909 = 1 in 11, ceiling 11.00x -- {0.0909/br:,.0f}x apart)")
order=np.argsort(-sc); rng=np.random.default_rng(SEED); N=len(sc)
print(f"\n  {'panel':>7} {'sites':>13} {'mutations found':>16} {'capture':>9} {'enrich':>8} {'RANDOM':>8}")
for k in (0.01,0.1,1.0,5.0,10.0):
    n=max(1,int(N*k/100)); c_=int(pos[order[:n]].sum())
    r_=int(pos[rng.permutation(N)[:n]].sum())
    print(f"  {k:6.2f}% {n:13,} {c_:16,} {100*c_/tot_mut:8.2f}% "
          f"{(100*c_/tot_mut)/k:8.3f} {(100*r_/tot_mut)/k:8.3f}")
