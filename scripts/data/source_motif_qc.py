"""Motif-QC: per-source deaminase signature vs expected; diagnose coordinate bugs."""
import numpy as np, pandas as pd
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
fa=Fasta(HG38)
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CHROMS=set([f'chr{i}' for i in range(1,23)]+['chrX','chrY'])

def trinuc(chrom,pos,ref,alt,base,d=0):
    """C- or A-centered trinuc at edited base shifted by offset d (strand-aware via ref/alt)."""
    p=pos+d
    try: w=fa[chrom][p-2:p+1].seq.upper()
    except: return None
    if len(w)!=3 or 'N' in w: return None
    if base=='C':
        if (ref,alt)==('C','T'): return w if w[1]=='C' else None
        if (ref,alt)==('G','A'):
            r=rc(w); return r if r[1]=='C' else None
    if base=='A':
        if (ref,alt)==('A','G'): return w if w[1]=='A' else None
        if (ref,alt)==('T','C'):
            r=rc(w); return r if r[1]=='A' else None
    return None

df=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
df['src']=df.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)
df=df[df.chrom.isin(CHROMS)]
df['cls']=df.apply(lambda r:'CBE' if (r.ref,r.alt) in [('C','T'),('G','A')] else ('ABE' if (r.ref,r.alt) in [('A','G'),('T','C')] else 'x'),axis=1)

EXPECT={'Doman_BE4_pilot':('CBE','rAPOBEC1->TpC'),'McGrath_2019_iPSC':('CBE','rAPOBEC1->TpC'),
        'Lei_2021':('CBE','rAPOBEC1(BE4max)->TpC'),'Yu_2020':('CBE','eng. YE1'),
        'Chen_2023_tdCBE':('CBE','TadA-CBE'),'Richter_2020':('ABE','TadA->TpA')}
def sig(sub,base,d=0):
    cs=[trinuc(r.chrom,int(r.pos),r.ref,r.alt,base,d) for r in sub.itertuples()]; cs=[c for c in cs if c]
    if not cs: return None
    s=pd.Series(cs).value_counts(normalize=True)
    key='T'+base
    yp=s[[k for k in s.index if k[0]=='T']].sum()
    return len(cs), yp, s.nlargest(3).to_dict()

print("=== as-is signature per source (d=0) ===")
for src,(cl,exp) in EXPECT.items():
    sub=df[(df.src==src)&(df.cls==cl)]; base='C' if cl=='CBE' else 'A'
    r=sig(sub,base)
    if r: 
        flag='  <-- FLAG' if (('TpC' in exp or 'TpA' in exp) and r[1]<0.4) else ''
        print(f"  {src:20s} {cl} exp={exp:24s} n={r[0]:5d}  T{base}={r[1]:.2f}  top={r[2]}{flag}")

print("\n=== Lei offset/strand scan (find the coordinate fix) ===")
lei=df[(df.src=='Lei_2021')&(df.cls=='CBE')]
for d in [-2,-1,0,1,2]:
    r=sig(lei,'C',d)
    if r: print(f"  offset d={d:+d}: n={r[0]:5d}  TpC={r[1]:.2f}  top={r[2]}")
# also try treating reported strand as flipped (swap ref/alt complement) at d=0
lei_flip=lei.copy()
lei_flip['ref']=lei['ref'].map(COMP); lei_flip['alt']=lei['alt'].map(COMP)
r=sig(lei_flip,'C',0)
if r: print(f"  strand-flip d=0: n={r[0]:5d}  TpC={r[1]:.2f}  top={r[2]}")
