#!/usr/bin/env python
"""THE SAME ATTACK THAT KILLED SECTION 31, NOW AIMED AT SECTION 35.

For the EDITOR I raced the trained model against `stem` alone -- one integer, no training, no
cancer data -- and stem WON at the panel sizes that matter. That retraction (section 32) is
the reason section 31's "cancer-trained model transfers" framing is gone.

I NEVER RAN THAT COMPARISON FOR THE CANCER GENOMIC RESULT. If stem alone also reaches ~6.8x
genomically in cancer, then section 35 says nothing about a trained model either: the honest
statement collapses to "hairpin geometry concentrates APOBEC mutations", the model is
decoration, and the deliverable is an annotation rather than a model.

Rank the whole genome by:
    the out-of-fold trained model   (needs PCAWG, needs training)
    stem alone                      (one integer from the universe file, NOTHING learned)
    hp_score alone                  (one float, nothing learned)
    random                          (the null)
and compare capture of the same 92,823 PCAWG mutations.
"""
import numpy as np, time
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
chroms=sorted(set(ctr.tolist()),key=lambda c:(len(c),c)); folds=[chroms[i::5] for i in range(5)]
s=np.load(f"{PC}/snvs.npz",allow_pickle=True)
keep=np.isin(s["donor"],list(train_donors))
mu_chrom=s["chrom"][keep]; mu_pos=s["pos"][keep].astype(np.int64)

mdl=[]; stem=[]; hp=[]; pos=[]; tot=0
for hold in folds:
    tr=~np.isin(ctr,hold)
    m=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.1,early_stopping=True,
                                     validation_fraction=0.1,random_state=SEED).fit(Xtr[tr],ytr[tr])
    for c in hold:
        U=np.load(f"{FEAT}/universe_chr{c}.npz")
        upos=U["pos1"] if "pos1" in U.files else U["pos"].astype(np.int64)+1
        X=np.column_stack([U[f] for f in FEATS]).astype(np.float32)
        mdl.append(m.predict_proba(X)[:,1].astype(np.float32))
        stem.append(U["stem"].astype(np.float32)); hp.append(U["hp_score"].astype(np.float32))
        mp=np.sort(mu_pos[mu_chrom==c]); i=np.clip(np.searchsorted(upos,mp),0,len(upos)-1)
        hit=np.zeros(len(upos),bool); hit[i[upos[i]==mp]]=True
        tot+=int(hit.sum()); pos.append(hit); del U,X
mdl=np.concatenate(mdl); stem=np.concatenate(stem); hp=np.concatenate(hp); pos=np.concatenate(pos)
N=len(pos); log(f"{N:,} sites, {tot:,} mutations")
rng=np.random.default_rng(SEED)
def cap(score,label):
    o=np.argsort(-(score+rng.uniform(0,1e-6,N).astype(np.float32)))
    out=[]
    for k in (0.01,0.1,1.0,5.0,10.0):
        n=max(1,int(N*k/100)); c_=int(pos[o[:n]].sum())
        out.append((k,c_,100*c_/tot,(100*c_/tot)/k))
    return out
res={n:cap(v,n) for v,n in ((mdl,"OOF model"),(stem,"stem ALONE"),(hp,"hp_score ALONE"))}
print(f"\n  {'ranking':14s} {'panel':>7} {'found':>9} {'capture':>9} {'enrich':>8}")
for name in ("OOF model","stem ALONE","hp_score ALONE"):
    for k,c_,cp,e in res[name]:
        print(f"  {name:14s} {k:6.2f}% {c_:9,} {cp:8.2f}% {e:8.3f}")
print(f"\n  {'panel':>7} {'model':>9} {'stem':>9} {'hp_score':>9} {'model - best untrained':>24}")
for i,k in enumerate((0.01,0.1,1.0,5.0,10.0)):
    em=res["OOF model"][i][3]; es=res["stem ALONE"][i][3]; eh=res["hp_score ALONE"][i][3]
    print(f"  {k:6.2f}% {em:9.3f} {es:9.3f} {eh:9.3f} {em-max(es,eh):24.3f}")
print("\n  If the model does not beat the untrained rankings, section 35 is about hairpin")
print("  geometry and not about a model -- exactly as section 32 concluded for the editor.")
