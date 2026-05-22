"""F1: full model vs simple probe. Feature ablation on the SAME held-out transfer test
(train Doman -> test McGrath, rAPOBEC1) + LOSO. Uses cached HyenaDNA(256) + handcrafted(81)
+ ENCODE g(4) + seq+-10. NOTE: union-pool negatives are Strategy-B (trinuc-matched) -> motif
partly removed; absolute AUROC differs from the random-neg 0.82. Comparison across feature sets
is fair (negatives fixed). torch-MPS, numpy AUROC."""
import json, numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu"); print("device",DEV)
RNG=np.random.default_rng(0); torch.manual_seed(0)
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
def auroc(y,s):
    y=np.asarray(y); s=np.asarray(s); o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(1,len(s)+1)
    ss=s[o]; i=0
    while i<len(s):
        j=i
        while j+1<len(s) and ss[j+1]==ss[i]: j+=1
        if j>i: r[o[i:j+1]]=(i+1+j+1)/2
        i=j+1
    p=y.sum(); n=len(y)-p
    return float('nan') if p==0 or n==0 else float((r[y==1].sum()-p*(p+1)/2)/(p*n))

df=pd.read_parquet('data/processed/train_table_v1.parquet')
cols=json.load(open('data/processed/feature_cols_v1.json')); FC,GC=cols['f'],cols['g']
df=df[df.cls=='CBE'].reset_index(drop=True)
# HyenaDNA
hy=torch.load('/tmp/hyenadna_union.pt',map_location='cpu')
def key(r): return f"{r['chrom']}:{int(r['pos'])}:{r['strand']}:{r['ref']}:{r['alt']}"
df['k']=df.apply(key,axis=1)
have=df.k.isin(hy); print(f"HyenaDNA coverage: {have.mean():.1%} ({have.sum()}/{len(df)})")
df=df[have].reset_index(drop=True)
HYE=np.stack([hy[k].numpy() for k in df.k]).astype(np.float32)
# seq +-10 one-hot (oriented to edited C)
def seq10(r):
    p=int(r['pos'])
    try: w=fa[r['chrom']][p-1-W:p+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if (r['ref'],r['alt'])==('G','A'): w=rc(w)
    return w
ss=[seq10(r) for _,r in df.iterrows()]
okmask=np.array([s is not None for s in ss]); df=df[okmask].reset_index(drop=True); HYE=HYE[okmask]
ss=[s for s in ss if s]
SEQ=np.zeros((len(ss),(2*W+1)*4),np.float32)
for i,s in enumerate(ss):
    for j,b in enumerate(s): SEQ[i,j*4+M[b]]=1
SEQ=SEQ[:,[c for c in range(SEQ.shape[1]) if c//4!=W]]
HAND=df[FC].values.astype(np.float32); GG=df[GC].values.astype(np.float32)
y=df.is_positive.values.astype(float); src=df.src.values
print(f"final CBE rows: {len(df)} (pos {int(y.sum())}/neg {int((1-y).sum())})")

FSETS={'seq10':SEQ,'hand':HAND,'hand+g':np.c_[HAND,GG],'hyena':HYE,'all':np.c_[SEQ,HAND,GG,HYE]}
class MLP(nn.Module):
    def __init__(s,d): super().__init__(); s.n=nn.Sequential(nn.Linear(d,128),nn.ReLU(),nn.Dropout(.3),nn.Linear(128,32),nn.ReLU(),nn.Dropout(.3),nn.Linear(32,1))
    def forward(s,x): return s.n(x).squeeze(-1)
def fit_eval(Xtr,ytr,Xte,yte,ep=80):
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6; Xtr=(Xtr-mu)/sd; Xte=(Xte-mu)/sd
    Xtr=torch.tensor(Xtr,device=DEV); yt=torch.tensor(ytr,dtype=torch.float32,device=DEV); Xte=torch.tensor(Xte,device=DEV)
    torch.manual_seed(0); m=MLP(Xtr.shape[1]).to(DEV); opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/max((ytr==1).sum(),1)],dtype=torch.float32,device=DEV); lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    for _ in range(ep): m.train(); opt.zero_grad(); lf(m(Xtr),yt).backward(); opt.step()
    m.eval()
    with torch.no_grad(): return auroc(yte,torch.sigmoid(m(Xte)).cpu().numpy())

negidx=np.where(y==0)[0]; RNG.shuffle(negidx); cut=int(0.2*len(negidx)); nte,ntr=negidx[:cut],negidx[cut:]
SRCS=['Doman_BE4_pilot','McGrath_2019_iPSC','Yu_2020','Chen_2023_tdCBE']
print("\n=== feature ablation: transfer AUROC (train->test), rAPOBEC1 pair + LOSO mean ===")
print(f"{'featureset':12s} {'Doman->McGr':>12s} {'McGr->Doman':>12s} {'LOSO mean':>10s}")
for name,X in FSETS.items():
    dm=fit_eval(X[np.r_[np.where(src=='McGrath_2019_iPSC')[0],ntr]],y[np.r_[np.where(src=='McGrath_2019_iPSC')[0],ntr]],
                X[np.r_[np.where(src=='Doman_BE4_pilot')[0],nte]],y[np.r_[np.where(src=='Doman_BE4_pilot')[0],nte]])  # train McGr->test Doman
    md=fit_eval(X[np.r_[np.where(src=='Doman_BE4_pilot')[0],ntr]],y[np.r_[np.where(src=='Doman_BE4_pilot')[0],ntr]],
                X[np.r_[np.where(src=='McGrath_2019_iPSC')[0],nte]],y[np.r_[np.where(src=='McGrath_2019_iPSC')[0],nte]])
    loso=[]
    for ho in SRCS:
        tr=np.r_[np.where((src!=ho)&(y==1))[0],ntr]; te=np.r_[np.where((src==ho)&(y==1))[0],nte]
        loso.append(fit_eval(X[tr],y[tr],X[te],y[te]))
    print(f"{name:12s} {md:12.3f} {dm:12.3f} {np.nanmean(loso):10.3f}")
print("\n[note] Strategy-B (trinuc-matched) negatives -> motif removed; for vs-0.82 (random negs) see F1b.")
