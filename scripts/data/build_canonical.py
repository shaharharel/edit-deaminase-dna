"""Canonical cleanup: v2 coords + genome-based strand resolution (recovers both strands) +
Lei -2 fix + per-source motif-QC. Output: canonical_sites.parquet (oriented to edited base)."""
import pandas as pd, numpy as np
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=set([f'chr{i}' for i in range(1,23)]+['chrX','chrY'])
v2=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
v2['src']=v2.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)
v2=v2[v2.chrom.isin(CH)].copy()
# Lei -2 coordinate fix
v2['pos']=v2.pos.astype(int)
v2.loc[v2.src=='Lei_2021','pos']-=2
def edit_class(r):
    p=(r.ref,r.alt)
    if p in [('C','T'),('G','A')]: return 'CBE'
    if p in [('A','G'),('T','C')]: return 'ABE'
    return None
fam_map={'Doman_BE4_pilot':'rAPOBEC1','McGrath_2019_iPSC':'rAPOBEC1','Lei_2021':'rAPOBEC1',
         'Yu_2020':'rAPOBEC1_eng','Chen_2023_tdCBE':'TadA_CBE','Richter_2020':'TadA','CHANGE-seq-BE':'mixed'}
rows=[]
for r in v2.itertuples():
    cl=edit_class(r)
    if cl is None: continue
    p=int(r.pos)
    try: w=fa[r.chrom][p-1-W:p+W].seq.upper()
    except: continue
    if len(w)!=2*W+1 or 'N' in w: continue
    center=w[W]
    if cl=='CBE':
        if center=='C': strand,ctx='+',w
        elif center=='G': strand,ctx='-',rc(w)
        else: continue
        eb='C'
    else:  # ABE
        if center=='A': strand,ctx='+',w
        elif center=='T': strand,ctx='-',rc(w)
        else: continue
        eb='A'
    rows.append((r.chrom,p,strand,cl,eb,ctx[W-1]+ctx[W]+ctx[W+1],ctx,r.src,fam_map.get(r.src,'other')))
df=pd.DataFrame(rows,columns=['chrom','pos','strand','edit_class','edit_base','trinuc','context','src','family'])
print(f"canonical positives: {len(df)} / {len(v2)} v2 ({len(df)/len(v2):.1%} oriented)")
print("by class:", df.edit_class.value_counts().to_dict())
df.to_parquet('data/processed/canonical_sites.parquet',index=False)
# motif-QC per source
print("\n=== motif-QC after cleanup (TpC for CBE / TpA for ABE) ===")
for src in df.src.unique():
    s=df[df.src==src]; cl=s.edit_class.mode()[0]; eb='C' if cl=='CBE' else 'A'
    tp=(s.trinuc.str[0]=='T').mean()
    top=s.trinuc.value_counts(normalize=True).head(3).to_dict()
    flag='' if tp>0.4 or s.family.iloc[0] in('rAPOBEC1_eng','TadA_CBE','mixed') else '  <-- LOW'
    print(f"  {src:22s} {cl} n={len(s):6d}  Tp{eb}={tp:.2f}  top={ {k:round(v,2) for k,v in top.items()} }{flag}")
