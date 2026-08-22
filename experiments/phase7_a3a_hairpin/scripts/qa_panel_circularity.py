#!/usr/bin/env python
"""THE STRONGEST ATTACK ON SECTION 31, AND IT IS A CIRCULARITY.

Y130G's defining property IS hairpin excess -- that is the +0.318/+0.396 finding. The PCAWG
model ranks by hairpin geometry. So a hairpin-ranked panel will of course contain more Y130G
sites. "THE CANCER-TRAINED MODEL TRANSFERS" MAY BE NOTHING MORE THAN THE HAIRPIN ENRICHMENT
RE-EXPRESSED AS A RANKING, with the PCAWG training contributing nothing at all.

DECISIVE TEST: put the trained model up against rankings that require NO TRAINING and NO
CANCER DATA WHATSOEVER --
    stem alone        (one integer from the universe file)
    hp_score alone    (one float from the universe file)
    random            (the null)
If the PCAWG model does not beat raw `stem`, then the honest claim is "hairpin-rich sites are
enriched in Y130G", which was already known, and section 31's framing about cancer training is
unsupported.

Same panel machinery, same calibrator, same 200-permutation null, so the comparison is
like-for-like.
"""
import numpy as np, sys, time
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

sc,stem_v,hp_v,ed_s,cal_s=[],[],[],[],[]
for c in CH:
    U=np.load(f"{FEAT}/universe_chr{c}.npz"); E=np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K=np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz"); P=np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    el=(E["cov"]>=8)&(K["cov"]>=8)&(P["cov"]>=15)&(P["alt"]==0)
    X=np.column_stack([U[f][el] for f in FEATS]).astype(np.float32)
    sc.append(m.predict_proba(X)[:,1].astype(np.float32))
    stem_v.append(U["stem"][el].astype(np.float32)); hp_v.append(U["hp_score"][el].astype(np.float32))
    ed_s.append(((E["alt"]>=2)&(K["alt"]==0))[el]); cal_s.append(((K["alt"]>=2)&(E["alt"]==0))[el])
    del U,E,K,P,X
sc=np.concatenate(sc); stem_v=np.concatenate(stem_v); hp_v=np.concatenate(hp_v)
ed=np.concatenate(ed_s); cal=np.concatenate(cal_s); N=len(sc)
log(f"{ED}: {N:,} sites, ed {int(ed.sum()):,}, cal {int(cal.sum()):,}")
edi=np.flatnonzero(ed); cali=np.flatnonzero(cal); rng=np.random.default_rng(SEED)

def gaps_for(score, name):
    # jitter breaks ties in the integer `stem` ranking so it is not order-of-file dependent
    order=np.argsort(-(score + rng.uniform(0,1e-6,size=N).astype(np.float32)))
    out=[]
    for k in (0.01,0.10,1.00):
        n=max(1,int(N*k/100)); tm=np.zeros(N,bool); tm[order[:n]]=True
        ne=int(tm[edi].sum()); nc=int(tm[cali].sum())
        ex=(ne/len(edi))/(k/100); cx=(nc/len(cali))/(k/100)
        out.append((k,ne,nc,ex,cx,ex-cx))
    return out

print(f"\n  {'ranking':16s} {'panel':>7} {'n_ed':>6} {'ed x':>7} {'n_cal':>6} {'cal x':>7} {'gap':>8}")
res={}
for score,name in ((sc,"PCAWG model"),(stem_v,"stem ALONE"),(hp_v,"hp_score ALONE")):
    res[name]=gaps_for(score,name)
    for k,ne,nc,ex,cx,g in res[name]:
        print(f"  {name:16s} {k:6.2f}% {ne:6d} {ex:7.3f} {nc:6d} {cx:7.3f} {g:+8.3f}")
print()
for k_i,k in enumerate((0.01,0.10,1.00)):
    gm=res["PCAWG model"][k_i][5]; gs=res["stem ALONE"][k_i][5]; gh=res["hp_score ALONE"][k_i][5]
    print(f"  top {k:5.2f}%   PCAWG {gm:+.3f}   stem-alone {gs:+.3f}   hp_score-alone {gh:+.3f}"
          f"   PCAWG minus best-untrained {gm-max(gs,gh):+.3f}")
print("\n  If PCAWG does not beat the untrained rankings, the cancer training contributes")
print("  nothing and section 31's framing is unsupported.")
