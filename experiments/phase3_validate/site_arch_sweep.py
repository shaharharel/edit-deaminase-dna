"""Autonomous iteration: architecture sweep + feature ablation + HARDER (trinuc-matched) background.
Within-editor held-out-chromosome. Arch: MLP / deep-MLP / CNN(local seq). Honest baselines."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=25; RNG=np.random.default_rng(0); DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
bn=pd.read_parquet('/tmp/bins_1mb_v3.parquet')
acc={(r.chrom,int(r.bin)):np.array([r.dnase,r.atac,r.rloop,r.h3k27ac],np.float32) for r in bn.itertuples()}
amu=np.nanmean(list(acc.values()),0); asd=np.nanstd(list(acc.values()),0)+1e-6
def ctxget(chrom,pos):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if w[W]=='G': w=rc(w)
    return w if w[W]=='C' else None
def trinuc(w): return w[W-1]+w[W]+w[W+1]
def onehot2d(w):
    a=np.zeros((4,2*W+1),np.float32)
    for j,b in enumerate(w):
        if b in M: a[M[b],j]=1
    return a
def accof(chrom,pos):
    a=acc.get((chrom,int(pos//1_000_000)))
    return (a-amu)/asd if (a is not None and not np.any(np.isnan(a))) else np.zeros(4,np.float32)
def site_recs(df):
    out=[]
    for r in df.itertuples():
        w=r.context if len(r.context)==2*10+1 else None  # canonical context is +-10; re-extract for +-25
        w=ctxget(r.chrom,int(r.pos))
        if w: out.append((r.chrom,int(r.pos),w,trinuc(w)))
    return out
def bg_matched(trinuc_counts, mult=4):
    """trinuc-matched random-C background: sample to match positives' trinuc distribution."""
    want={t:c*mult for t,c in trinuc_counts.items()}; got={t:0 for t in want}; out=[]
    wt=np.array([CL[c] for c in CH],float); wt/=wt.sum(); t=0; need=sum(want.values())
    while len(out)<need and t<need*60:
        t+=1;c=RNG.choice(CH,p=wt);p=int(RNG.integers(W+2,CL[c]-W-2));w=ctxget(c,p)
        if not w: continue
        tn=trinuc(w)
        if got.get(tn,0)<want.get(tn,0): out.append((c,p,w,tn)); got[tn]=got.get(tn,0)+1
    return out
def auroc(y,sc):
    o=np.argsort(sc);r=np.empty(len(sc));r[o]=np.arange(1,len(sc)+1);p=y.sum();n=len(y)-p
    return float((r[y==1].sum()-p*(p+1)/2)/(p*n))
def rec(y,sc,K=0.10):
    o=np.argsort(-sc);return float(y[o[:int(K*len(y))]].sum()/y.sum())
class MLPm(nn.Module):
    def __init__(s,d,deep=False):
        super().__init__()
        if deep: s.n=nn.Sequential(nn.Linear(d,128),nn.ReLU(),nn.Dropout(.3),nn.Linear(128,64),nn.ReLU(),nn.Dropout(.3),nn.Linear(64,32),nn.ReLU(),nn.Linear(32,1))
        else: s.n=nn.Sequential(nn.Linear(d,64),nn.ReLU(),nn.Dropout(.3),nn.Linear(64,1))
    def forward(s,x): return s.n(x).squeeze(-1)
class CNN(nn.Module):
    def __init__(s,acc_d=4):
        super().__init__()
        s.c=nn.Sequential(nn.Conv1d(4,32,7,padding=3),nn.ReLU(),nn.Conv1d(32,32,5,padding=2),nn.ReLU(),nn.AdaptiveMaxPool1d(1))
        s.h=nn.Sequential(nn.Linear(32+acc_d,32),nn.ReLU(),nn.Dropout(.3),nn.Linear(32,1))
    def forward(s,seq,a): return s.h(torch.cat([s.c(seq).squeeze(-1),a],1)).squeeze(-1)
def fit(Xtr,ytr,Xte,deep=False,ep=40):
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6;Xtr=(Xtr-mu)/sd;Xte=(Xte-mu)/sd
    m=MLPm(Xtr.shape[1],deep).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/max((ytr==1).sum(),1)],dtype=torch.float32,device=DEV);lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(ytr,dtype=torch.float32,device=DEV);idx=np.arange(len(Xs));torch.manual_seed(0)
    for e in range(ep):
        RNG.shuffle(idx)
        for i in range(0,len(idx),4096):
            bi=idx[i:i+4096];m.train();opt.zero_grad();lf(m(Xs[bi]),ys[bi]).backward();opt.step()
    m.eval()
    with torch.no_grad(): return torch.sigmoid(m(torch.tensor(Xte,device=DEV))).cpu().numpy()
s=pd.read_parquet('data/processed/canonical_sites.parquet');s=s[s.chrom.isin(CH)]
trainchr=set(CH[::2])
for ed in ['Lei_2021','Yu_2020']:
    P=site_recs(s[(s.edit_class=='CBE')&(s.src==ed)])
    if len(P)<300: continue
    from collections import Counter
    tc=Counter(p[3] for p in P)
    B=bg_matched(tc,4)  # trinuc-MATCHED background (harder)
    rows=P+B; y=np.r_[np.ones(len(P)),np.zeros(len(B))]
    chrom=np.array([r[0] for r in rows]); tr=np.isin(chrom,list(trainchr))
    seq=np.stack([onehot2d(r[2]) for r in rows]); a=np.stack([accof(r[0],r[1]) for r in rows])
    mot10=seq[:, :, W-10:W+11].reshape(len(rows),-1)  # +-10 flatten for MLP
    print(f"\n=== {ed}: within-editor held-out-chrom, TRINUC-MATCHED bg (pos {len(P)}, bg {len(B)}) ===")
    # feature ablations (MLP)
    for name,X in [('motif(+-10)',mot10),('acc',a),('motif+acc',np.c_[mot10,a])]:
        sc=fit(X[tr],y[tr],X[~tr]); print(f"  MLP {name:14s}  AUROC {auroc(y[~tr],sc):.3f}  R@10% {rec(y[~tr],sc):.1%}")
    sc=fit(np.c_[mot10,a][tr],y[tr],np.c_[mot10,a][~tr],deep=True); print(f"  deepMLP motif+acc  AUROC {auroc(y[~tr],sc):.3f}  R@10% {rec(y[~tr],sc):.1%}")
    # CNN over +-25 seq + acc
    Strn=torch.tensor(seq[tr],device=DEV);atrn=torch.tensor((a[tr]),device=DEV);ytr=y[tr]
    m=CNN().to(DEV);opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/max((ytr==1).sum(),1)],dtype=torch.float32,device=DEV);lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    ys=torch.tensor(ytr,dtype=torch.float32,device=DEV);idx=np.arange(len(Strn));torch.manual_seed(0)
    for e in range(40):
        RNG.shuffle(idx)
        for i in range(0,len(idx),2048):
            bi=idx[i:i+2048];m.train();opt.zero_grad();lf(m(Strn[bi],atrn[bi]),ys[bi]).backward();opt.step()
    m.eval()
    with torch.no_grad(): sc=torch.sigmoid(m(torch.tensor(seq[~tr],device=DEV),torch.tensor(a[~tr],device=DEV))).cpu().numpy()
    print(f"  CNN(+-25)+acc      AUROC {auroc(y[~tr],sc):.3f}  R@10% {rec(y[~tr],sc):.1%}")
