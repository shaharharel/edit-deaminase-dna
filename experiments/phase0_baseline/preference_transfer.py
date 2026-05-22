"""P0: Is guide-independent deaminase preference learnable & transferable at site level?
Sequence-only local (+-W) context. Numpy-only (env scipy/sklearn broken)."""
import numpy as np, pandas as pd
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
W=10; RNG=np.random.default_rng(0)
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
fa=Fasta(HG38); CHROMS=[f'chr{i}' for i in range(1,23)]+['chrX']
CHROMLEN={c:len(fa[c]) for c in CHROMS}

def auroc(y,s):
    y=np.asarray(y); s=np.asarray(s); order=np.argsort(s); r=np.empty(len(s)); r[order]=np.arange(1,len(s)+1)
    # average ties
    s_sorted=s[order]; i=0
    while i<len(s):
        j=i
        while j+1<len(s) and s_sorted[j+1]==s_sorted[i]: j+=1
        if j>i:
            avg=(i+1+j+1)/2.0
            r[order[i:j+1]]=avg
        i=j+1
    np_=y.sum(); nn=len(y)-np_
    if np_==0 or nn==0: return float('nan')
    return (r[y==1].sum()-np_*(np_+1)/2)/(np_*nn)

def fit_lr(X,y,iters=400,lr=0.5,l2=1e-3):
    n,d=X.shape; w=np.zeros(d); b=0.0
    for _ in range(iters):
        z=X@w+b; p=1/(1+np.exp(-z)); g=p-y
        w-=lr*(X.T@g/n + l2*w); b-=lr*g.mean()
    return w,b
def proba(Xb,X): w,b=Xb; return 1/(1+np.exp(-(X@w+b)))

def ctx(chrom,pos,ref,alt,base):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if base=='C':
        if (ref,alt)==('C','T') and w[W]=='C': return w
        if (ref,alt)==('G','A'):
            r=rc(w); return r if r[W]=='C' else None
    if base=='A':
        if (ref,alt)==('A','G') and w[W]=='A': return w
        if (ref,alt)==('T','C'):
            r=rc(w); return r if r[W]=='A' else None
    return None
def sample_neg(n,base):
    out=[]; wt=np.array([CHROMLEN[c] for c in CHROMS],float); wt/=wt.sum(); tries=0
    while len(out)<n and tries<n*60:
        tries+=1; c=RNG.choice(CHROMS,p=wt); p=int(RNG.integers(W+2,CHROMLEN[c]-W-2))
        w=fa[c][p-1-W:p+W].seq.upper()
        if len(w)!=2*W+1 or 'N' in w: continue
        if w[W]==base: out.append(w)
        elif w[W]==COMP[base]: out.append(rc(w))
    return out
M={'A':0,'C':1,'G':2,'T':3}
def onehot(cs):
    X=np.zeros((len(cs),(2*W+1)*4),np.float32)
    for i,s in enumerate(cs):
        for j,b in enumerate(s): X[i,j*4+M[b]]=1
    keep=[k for k in range(X.shape[1]) if k//4 != W]  # drop center
    return X[:,keep]

df=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
df['src']=df.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True); df=df[df.chrom.isin(CHROMS)]
CBE=['Doman_BE4_pilot','McGrath_2019_iPSC','Yu_2020','Lei_2021','Chen_2023_tdCBE']
pos={}
for s in CBE:
    p=[ctx(r.chrom,int(r.pos),r.ref,r.alt,'C') for r in df[df.src==s].itertuples()]; p=[c for c in p if c]
    if len(p)>=300: pos[s]=p
print("CBE oriented-C positives:", {k:len(v) for k,v in pos.items()})
data={s:(onehot(p+sample_neg(len(p),'C')), np.r_[np.ones(len(p)),np.zeros(len(p))]) for s,p in pos.items()}
srcs=list(data.keys()); models={s:fit_lr(*data[s]) for s in srcs}
print("\n=== Transfer AUROC (rows=train, cols=test) ===")
print("train\\test    "+"  ".join(f"{s[:10]:>10s}" for s in srcs))
for tr in srcs:
    print(f"{tr[:12]:12s}  "+"  ".join(f"{auroc(data[te][1],proba(models[tr],data[te][0])):10.3f}" for te in srcs))
print("\n=== within-source 5-fold CV (ceiling) ===")
for s in srcs:
    X,y=data[s]; idx=RNG.permutation(len(y)); folds=np.array_split(idx,5); aucs=[]
    for f in folds:
        te=np.zeros(len(y),bool); te[f]=True
        m=fit_lr(X[~te],y[~te]); aucs.append(auroc(y[te],proba(m,X[te])))
    print(f"  {s[:14]:14s} {np.mean(aucs):.3f} +- {np.std(aucs):.3f}")
if 'Doman_BE4_pilot' in data and 'McGrath_2019_iPSC' in data:
    Xa,ya=data['Doman_BE4_pilot']; Xb,yb=data['McGrath_2019_iPSC']; m=fit_lr(np.r_[Xa,Xb],np.r_[ya,yb])
    print("\n=== canonical rAPOBEC1 (Doman+McGrath) -> test ===")
    for te in srcs: print(f"  -> {te[:14]:14s} {auroc(data[te][1],proba(m,data[te][0])):.3f}")
