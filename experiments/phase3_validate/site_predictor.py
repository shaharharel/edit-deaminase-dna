"""SITE-LEVEL off-target predictor. Each candidate C: local motif (+-10 one-hot, center) + regional
accessibility (1Mb bin tracks). Train on Lei (CBE), recover HELD-OUT Yu sites among random-C background.
Honest baselines: gene-density / accessibility / motif / full."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; RNG=np.random.default_rng(0); DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
bn=pd.read_parquet('/tmp/bins_1mb_v3.parquet')
acc={(r.chrom,int(r.bin)):np.array([r.dnase,r.atac,r.rloop,r.h3k27ac],np.float32) for r in bn.itertuples()}
accarr=np.array([v for v in acc.values()]); amu=np.nanmean(accarr,0); asd=np.nanstd(accarr,0)+1e-6
rgd=pd.read_csv('/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt',sep='\t',header=None,usecols=[2,4,5,12],names=['chrom','ts','te','n2'])
rgd['bb']=((rgd.ts+rgd.te)//2)//1_000_000; GD=rgd.groupby(['chrom','bb']).n2.nunique().to_dict()
def onehot(w):
    oh=np.zeros((2*W+1)*4,np.float32)
    for j,b in enumerate(w):
        if b in M: oh[j*4+M[b]]=1
    return oh[[k for k in range(len(oh)) if k//4!=W]]
def site_row(chrom,pos,context):
    if len(context)!=2*W+1 or 'N' in context or context[W]!='C': return None
    b1=int(pos//1_000_000); a=acc.get((chrom,b1))
    if a is None or np.any(np.isnan(a)): a=amu
    return onehot(context),a,np.float32(GD.get((chrom,b1),0))
def from_sites(df):
    O=[];A=[];G=[]
    for r in df.itertuples():
        x=site_row(r.chrom,int(r.pos),r.context)
        if x: O.append(x[0]);A.append(x[1]);G.append(x[2])
    return np.array(O),np.array(A),np.array(G)
def rand_C(n):
    O=[];A=[];G=[];wt=np.array([CL[c] for c in CH],float);wt/=wt.sum();t=0
    while len(O)<n and t<n*40:
        t+=1;c=RNG.choice(CH,p=wt);p=int(RNG.integers(W+2,CL[c]-W-2))
        try: w=fa[c][p-1-W:p+W].seq.upper()
        except: continue
        if len(w)!=2*W+1 or 'N' in w: continue
        if w[W]=='G': w=rc(w)
        if w[W]!='C': continue
        x=site_row(c,p,w)
        if x: O.append(x[0]);A.append(x[1]);G.append(x[2])
    return np.array(O),np.array(A),np.array(G)
s=pd.read_parquet('data/processed/canonical_sites.parquet'); s=s[s.chrom.isin(CH)]
lei=s[(s.edit_class=='CBE')&(s.src=='Lei_2021')]; yu=s[(s.edit_class=='CBE')&(s.src=='Yu_2020')]
print("extracting features...")
Lo,La,Lg=from_sites(lei); Yo,Ya,Yg=from_sites(yu); Bo,Ba,Bg=rand_C(40000); To,Ta,Tg=rand_C(200000)
print(f"train: Lei {len(Lo)} pos + {len(Bo)} bg | test: Yu {len(Yo)} pos + {len(To)} bg")
def auroc(y,sc):
    o=np.argsort(sc);r=np.empty(len(sc));r[o]=np.arange(1,len(sc)+1);p=y.sum();n=len(y)-p
    return float((r[y==1].sum()-p*(p+1)/2)/(p*n))
def recall_at(y,sc,Ks=(0.001,0.01,0.05,0.10)):
    o=np.argsort(-sc);tot=y.sum();return {K:float(y[o[:int(K*len(y))]].sum()/tot) for K in Ks}
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,64),nn.ReLU(),nn.Dropout(.3),nn.Linear(64,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def build(O,A,G,which):
    c=[]
    if 'motif' in which: c.append(O)
    if 'acc' in which: c.append((A-amu)/asd)
    if 'gd' in which: c.append(G.reshape(-1,1))
    return np.c_[tuple(c)].astype(np.float32)
def run(which):
    Xtr=build(np.r_[Lo,Bo],np.r_[La,Ba],np.r_[Lg,Bg],which); ytr=np.r_[np.ones(len(Lo)),np.zeros(len(Bo))]
    Xte=build(np.r_[Yo,To],np.r_[Ya,Ta],np.r_[Yg,Tg],which); yte=np.r_[np.ones(len(Yo)),np.zeros(len(To))]
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6;Xtr=(Xtr-mu)/sd;Xte=(Xte-mu)/sd
    torch.manual_seed(0);m=MLP(Xtr.shape[1]).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/(ytr==1).sum()],dtype=torch.float32,device=DEV);lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(ytr,dtype=torch.float32,device=DEV);idx=np.arange(len(Xs))
    for ep in range(30):
        RNG.shuffle(idx)
        for i in range(0,len(idx),4096):
            bi=idx[i:i+4096];m.train();opt.zero_grad();lf(m(Xs[bi]),ys[bi]).backward();opt.step()
    m.eval()
    with torch.no_grad(): sc=torch.sigmoid(m(torch.tensor(Xte,device=DEV))).cpu().numpy()
    return auroc(yte,sc),recall_at(yte,sc)
print("\n=== SITE-LEVEL: held-out Yu (CBE) recovery vs random-C background (trained on Lei) ===")
print(f"{'features':16s} {'AUROC':>6s} {'R@0.1%':>7s} {'R@1%':>6s} {'R@5%':>6s} {'R@10%':>6s}")
for which in [['gd'],['acc'],['motif'],['motif','acc'],['motif','acc','gd']]:
    au,rcl=run(which)
    print(f"{'+'.join(which):16s} {au:6.3f} {rcl[0.001]:7.1%} {rcl[0.01]:6.1%} {rcl[0.05]:6.1%} {rcl[0.10]:6.1%}")
