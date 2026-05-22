"""P3: predicted INTRINSIC (sequence-only, f) per-gene CBE burden, to correlate vs BE4 gold index.
Train rAPOBEC1 preference model (Doman+McGrath edited C vs random genomic C), then score every
gene's CDS C's and aggregate mean predicted rate = predicted intrinsic susceptibility."""
import numpy as np, pandas as pd
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
REFGENE="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt"
W=10; RNG=np.random.default_rng(0); COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
fa=Fasta(HG38); CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
M={'A':0,'C':1,'G':2,'T':3}
def onehot(cs):
    X=np.zeros((len(cs),(2*W+1)*4),np.float32)
    for i,s in enumerate(cs):
        for j,b in enumerate(s): X[i,j*4+M[b]]=1
    return X[:,[k for k in range(X.shape[1]) if k//4!=W]]
def fit_lr(X,y,it=400,lr=0.5,l2=1e-3):
    n,d=X.shape; w=np.zeros(d); b=0.0
    for _ in range(it):
        p=1/(1+np.exp(-(X@w+b))); g=p-y; w-=lr*(X.T@g/n+l2*w); b-=lr*g.mean()
    return w,b
def proba(M_,X): w,b=M_; return 1/(1+np.exp(-(X@w+b)))
def ctx_edit(chrom,pos,ref,alt):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if (ref,alt)==('C','T'): return w if w[W]=='C' else None
    if (ref,alt)==('G','A'):
        r=rc(w); return r if r[W]=='C' else None
    return None
def sample_neg(n):
    out=[]; wt=np.array([CL[c] for c in CH],float); wt/=wt.sum(); t=0
    while len(out)<n and t<n*60:
        t+=1; c=RNG.choice(CH,p=wt); p=int(RNG.integers(W+2,CL[c]-W-2)); w=fa[c][p-1-W:p+W].seq.upper()
        if len(w)!=2*W+1 or 'N' in w: continue
        if w[W]=='C': out.append(w)
        elif w[W]=='G': out.append(rc(w))
    return out

# --- train rAPOBEC1 preference ---
df=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
df['src']=df.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)
rap=df[(df.src.isin(['Doman_BE4_pilot','McGrath_2019_iPSC']))]
pos=[ctx_edit(r.chrom,int(r.pos),r.ref,r.alt) for r in rap.itertuples()]; pos=[c for c in pos if c]
neg=sample_neg(len(pos))
X=onehot(pos+neg); y=np.r_[np.ones(len(pos)),np.zeros(len(neg))]
model=fit_lr(X,y); print(f"trained rAPOBEC1 preference on {len(pos)} pos / {len(neg)} neg")

# --- CDS C positions per gene (longest isoform) ---
rg=pd.read_csv(REFGENE,sep='\t',header=None,names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','name2','s1','s2','frames'])
rg=rg[rg.chrom.isin(CH)]; rg=rg[rg.cdsStart<rg.cdsEnd].copy(); rg['span']=rg.txEnd-rg.txStart
longest=rg.sort_values('span').groupby('name2').tail(1)
K=80; rows=[]
for r in longest.itertuples():
    ss=[int(x) for x in r.exonStarts.rstrip(',').split(',') if x]
    ee=[int(x) for x in r.exonEnds.rstrip(',').split(',') if x]
    cds=[]
    for a,b in zip(ss,ee):
        a2,b2=max(a,r.cdsStart),min(b,r.cdsEnd)
        if b2>a2: cds.append((a2,b2))
    # collect coding-strand C positions
    cpos=[]
    for a,b in cds:
        try: seq=fa[r.chrom][a:b].seq.upper()
        except: continue
        tb='C' if r.strand=='+' else 'G'
        for i,base in enumerate(seq):
            if base==tb: cpos.append(a+i+1)  # 1-based
    if len(cpos)<5: continue
    if len(cpos)>K: cpos=list(RNG.choice(cpos,K,replace=False))
    ctxs=[]
    for p in cpos:
        w=fa[r.chrom][p-1-W:p+W].seq.upper()
        if len(w)!=2*W+1 or 'N' in w: continue
        ctxs.append(w if r.strand=='+' else rc(w))
    if len(ctxs)<5: continue
    pr=proba(model,onehot(ctxs))
    gc=np.mean([ (s.count('G')+s.count('C'))/len(s) for s in ctxs ])  # GC control
    rows.append((r.name2,len(ctxs),float(pr.mean()),float(gc)))
out=pd.DataFrame(rows,columns=['gene','n_C_scored','pred_burden','gc_mean'])
out.to_parquet('data/processed/pred_gene_burden_v1.parquet',index=False)
print(f"scored {len(out)} genes -> data/processed/pred_gene_burden_v1.parquet")
print(out.sort_values('pred_burden',ascending=False).head(10).to_string(index=False))
