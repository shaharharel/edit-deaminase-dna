"""F1b: build a clean, current comparison dataset for full-model-vs-probe.
Positives = Doman+McGrath rAPOBEC1 CBE; negatives = random genomic C. Both C-oriented,
513bp context (for HyenaDNA) + the +-10 window derivable from it. Random-neg regime = 0.82-comparable."""
import numpy as np, pandas as pd
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=256; RNG=np.random.default_rng(0); COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
def ctx(chrom,pos,flip):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or w.count('N')>5: return None
    return rc(w) if flip else w

df=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
df['src']=df.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)
rap=df[df.src.isin(['Doman_BE4_pilot','McGrath_2019_iPSC']) & df.ref.isin(['C','G'])]
rows=[]
for r in rap.itertuples():
    flip = (r.ref,r.alt)==('G','A')
    s=ctx(r.chrom,int(r.pos),flip)
    if s and s[W]=='C': rows.append((f"{r.chrom}:{int(r.pos)}:{r.strand}:{r.ref}:{r.alt}",r.src,1,s))
npos=len(rows); print(f"positives (rAPOBEC1, C-oriented): {npos}")
# random genomic C negatives
neg=0; wt=np.array([CL[c] for c in CH],float); wt/=wt.sum(); tries=0
posset=set(k for k,_,_,_ in rows)
while neg<npos and tries<npos*80:
    tries+=1; c=RNG.choice(CH,p=wt); p=int(RNG.integers(W+2,CL[c]-W-2))
    b=fa[c][p-1:p].seq.upper()
    if b not in('C','G'): continue
    flip=(b=='G'); s=ctx(c,p,flip)
    if not s or s[W]!='C': continue
    k=f"{c}:{p}:{'+' if not flip else '-'}:C:T"
    if k in posset: continue
    rows.append((k,'random_neg',0,s)); neg+=1
print(f"negatives (random genomic C): {neg}")
out=pd.DataFrame(rows,columns=['key','src','label','seq'])
out.to_parquet('/tmp/compare_sites.parquet',index=False)
print(f"wrote /tmp/compare_sites.parquet  total {len(out)}  (pos {out.label.sum()}/neg {(out.label==0).sum()})")
print("seq len check:", out.seq.str.len().min(), out.seq.str.len().max())
