"""Within-editor, held-out-CHROMOSOME site-level recall (correct motif for the editor + accessibility).
The realistic 'predict THIS editor's off-targets in unseen genome' test. Yu and Lei separately."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; RNG=np.random.default_rng(0); DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
bn=pd.read_parquet('/tmp/bins_1mb_v3.parquet')
acc={(r.chrom,int(r.bin)):np.array([r.dnase,r.atac,r.rloop,r.h3k27ac],np.float32) for r in bn.itertuples()}
amu=np.nanmean(list(acc.values()),0); asd=np.nanstd(list(acc.values()),0)+1e-6
def oh(w):
    a=np.zeros((2*W+1)*4,np.float32)
    for j,b in enumerate(w):
        if b in M: a[j*4+M[b]]=1
    return a[[k for k in range(len(a)) if k//4!=W]]
def row(chrom,pos,ctx):
    if len(ctx)!=2*W+1 or 'N' in ctx or ctx[W]!='C': return None
    a=acc.get((chrom,int(pos//1_000_000)))
    if a is None or np.any(np.isnan(a)): a=amu
    return oh(ctx),(a-amu)/asd,chrom
def sites(df):
    O=[];A=[];C=[]
    for r in df.itertuples():
        x=row(r.chrom,int(r.pos),r.context)
        if x: O.append(x[0]);A.append(x[1]);C.append(x[2])
    return np.array(O),np.array(A),np.array(C)
def bg(n):
    O=[];A=[];C=[];wt=np.array([CL[c] for c in CH],float);wt/=wt.sum();t=0
    while len(O)<n and t<n*40:
        t+=1;c=RNG.choice(CH,p=wt);p=int(RNG.integers(W+2,CL[c]-W-2))
        try:w=fa[c][p-1-W:p+W].seq.upper()
        except:continue
        if len(w)!=2*W+1 or 'N' in w: continue
        if w[W]=='G': w=rc(w)
        if w[W]!='C': continue
        x=row(c,p,w)
        if x: O.append(x[0]);A.append(x[1]);C.append(x[2])
    return np.array(O),np.array(A),np.array(C)
def auroc(y,sc):
    o=np.argsort(sc);r=np.empty(len(sc));r[o]=np.arange(1,len(sc)+1);p=y.sum();n=len(y)-p
    return float((r[y==1].sum()-p*(p+1)/2)/(p*n))
def rec(y,sc,Ks=(0.001,0.01,0.05,0.10)):
    o=np.argsort(-sc);tot=y.sum();return {K:float(y[o[:int(K*len(y))]].sum()/tot) for K in Ks}
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,64),nn.ReLU(),nn.Dropout(.3),nn.Linear(64,1))
    def forward(s,x):return s.n(x).squeeze(-1)
s=pd.read_parquet('data/processed/canonical_sites.parquet');s=s[s.chrom.isin(CH)]
print("extracting background (shared)..."); Bo,Ba,Bc=bg(250000)
trainchr=set(CH[::2]); 
for ed in ['Yu_2020','Lei_2021','McGrath_2019_iPSC']:
    d=s[(s.edit_class=='CBE')&(s.src==ed)]
    Po,Pa,Pc=sites(d)
    if len(Po)<200: print(ed,'too few'); continue
    trm=np.isin(Pc,list(trainchr)); btr=np.isin(Bc,list(trainchr))
    Xtr=np.c_[np.r_[Po[trm],Bo[btr]], np.r_[Pa[trm],Ba[btr]]].astype(np.float32); ytr=np.r_[np.ones(trm.sum()),np.zeros(btr.sum())]
    Xte=np.c_[np.r_[Po[~trm],Bo[~btr]], np.r_[Pa[~trm],Ba[~btr]]].astype(np.float32); yte=np.r_[np.ones((~trm).sum()),np.zeros((~btr).sum())]
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6;Xtr=(Xtr-mu)/sd;Xte=(Xte-mu)/sd
    torch.manual_seed(0);m=MLP(Xtr.shape[1]).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/(ytr==1).sum()],dtype=torch.float32,device=DEV);lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(ytr,dtype=torch.float32,device=DEV);idx=np.arange(len(Xs))
    for ep in range(40):
        RNG.shuffle(idx)
        for i in range(0,len(idx),4096):
            bi=idx[i:i+4096];m.train();opt.zero_grad();lf(m(Xs[bi]),ys[bi]).backward();opt.step()
    m.eval()
    with torch.no_grad(): sc=torch.sigmoid(m(torch.tensor(Xte,device=DEV))).cpu().numpy()
    r=rec(yte,sc)
    print(f"  {ed:18s} (motif+acc, within-editor held-out-chrom): AUROC {auroc(yte,sc):.3f} | R@1% {r[0.01]:.1%} R@5% {r[0.05]:.1%} R@10% {r[0.10]:.1%}  (test pos {int(yte.sum())})")
