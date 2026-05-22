"""P2 v1: factorized f x g rate model. CBE leave-one-source-out + ABE CV.
Ablation: f-only vs f+g. Tests cross-source transfer + whether accessibility (g) adds.
torch-MPS, numpy AUROC (env scipy/sklearn broken)."""
import json, numpy as np, pandas as pd, torch, torch.nn as nn

DEV = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print("device:", DEV)
RNG = np.random.default_rng(0); torch.manual_seed(0)

df = pd.read_parquet("data/processed/train_table_v1.parquet")
cols = json.load(open("data/processed/feature_cols_v1.json"))
FCOLS, GCOLS = cols["f"], cols["g"]

def auroc(y, s):
    y=np.asarray(y); s=np.asarray(s); o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(1,len(s)+1)
    ss=s[o]; i=0
    while i<len(s):
        j=i
        while j+1<len(s) and ss[j+1]==ss[i]: j+=1
        if j>i: r[o[i:j+1]]=(i+1+j+1)/2.0
        i=j+1
    p=y.sum(); n=len(y)-p
    return float('nan') if p==0 or n==0 else float((r[y==1].sum()-p*(p+1)/2)/(p*n))

class MLP(nn.Module):
    def __init__(s,d):
        super().__init__(); s.net=nn.Sequential(nn.Linear(d,64),nn.ReLU(),nn.Dropout(0.3),
            nn.Linear(64,32),nn.ReLU(),nn.Dropout(0.3),nn.Linear(32,1))
    def forward(s,x): return s.net(x).squeeze(-1)

def fit_eval(Xtr,ytr,Xte,yte,epochs=60):
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6
    Xtr=(Xtr-mu)/sd; Xte=(Xte-mu)/sd
    Xtr=torch.tensor(Xtr,dtype=torch.float32,device=DEV); ytr_t=torch.tensor(ytr,dtype=torch.float32,device=DEV)
    Xte=torch.tensor(Xte,dtype=torch.float32,device=DEV)
    m=MLP(Xtr.shape[1]).to(DEV); opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/max((ytr==1).sum(),1)],dtype=torch.float32,device=DEV)
    lossf=nn.BCEWithLogitsLoss(pos_weight=pw)
    for _ in range(epochs):
        m.train(); opt.zero_grad(); out=m(Xtr); loss=lossf(out,ytr_t); loss.backward(); opt.step()
    m.eval()
    with torch.no_grad(): s=torch.sigmoid(m(Xte)).cpu().numpy()
    return auroc(yte,s)

def matrix(cls):
    sub=df[df.cls==cls].reset_index(drop=True)
    pos=sub[sub.is_positive]; neg=sub[~sub.is_positive]
    srcs=sorted([s for s,n in pos.src.value_counts().items() if n>=300])
    print(f"\n===== {cls} =====  pos sources={[(s,int((pos.src==s).sum())) for s in srcs]}  neg={len(neg)}")
    negidx=RNG.permutation(len(neg)); cut=int(0.2*len(neg))
    neg_te=neg.iloc[negidx[:cut]]; neg_tr=neg.iloc[negidx[cut:]]
    for use_g in [False,True]:
        feats=FCOLS+(GCOLS if use_g else [])
        tag="f+g" if use_g else "f  "
        if len(srcs)>=2:  # LOSO
            res=[]
            for ho in srcs:
                tr=pd.concat([pos[pos.src!=ho],neg_tr]); te=pd.concat([pos[pos.src==ho],neg_te])
                a=fit_eval(tr[feats].values,tr.is_positive.values.astype(float),
                           te[feats].values,te.is_positive.values.astype(float))
                res.append((ho,a))
            print(f"  [{tag}] LOSO AUROC: "+"  ".join(f"{h[:12]}={a:.3f}" for h,a in res)+f"   mean={np.mean([a for _,a in res]):.3f}")
        else:  # CV (single source)
            p=pos.sample(frac=1,random_state=0); folds=np.array_split(np.arange(len(p)),5)
            nf=np.array_split(RNG.permutation(len(neg)),5); aucs=[]
            for k in range(5):
                te=pd.concat([p.iloc[folds[k]],neg.iloc[nf[k]]])
                tr=pd.concat([p.drop(p.index[folds[k]]),neg.drop(neg.index[nf[k]])])
                aucs.append(fit_eval(tr[feats].values,tr.is_positive.values.astype(float),
                                     te[feats].values,te.is_positive.values.astype(float)))
            print(f"  [{tag}] 5-fold CV AUROC: {np.mean(aucs):.3f} +- {np.std(aucs):.3f}")

for cls in ["CBE","ABE"]:
    matrix(cls)
print("\n[done] v1 baseline (Strategy-B negatives). Next: 3 negative regimes + gold validation.")
