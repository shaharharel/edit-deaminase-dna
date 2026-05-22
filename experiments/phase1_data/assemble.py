"""P1: assemble local training table = labels + f (handcrafted) + g (encode), guide-indep only.
Drops Lei (coord bug) + CHANGE-seq (guide-dependent). Builds class (CBE/ABE) + regime negatives."""
import numpy as np, pandas as pd
KEY=['chrom','pos','strand']
f=pd.read_parquet('data/processed/handcrafted_features_v3.parquet')
g=pd.read_parquet('/tmp/encode_features_union.parquet')
lab=pd.read_parquet('data/processed/dna_offtarget_combined_B.parquet')
print('shapes f',f.shape,'g',g.shape,'lab',lab.shape)

def cls(r):
    p=(r['ref'],r['alt'])
    if p in [('C','T'),('G','A')]: return 'CBE'
    if p in [('A','G'),('T','C')]: return 'ABE'
    return 'x'
lab['cls']=lab.apply(cls,axis=1)
lab['src']=lab.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)

DROP_POS={'Lei_2021','CHANGE-seq-BE'}  # Lei coord bug; CHANGE-seq guide-dependent
# keep: all negatives (is_positive False), plus positives not in DROP_POS
keep = (~lab.is_positive) | (~lab.src.isin(DROP_POS))
lab=lab[keep & (lab.cls!='x')].copy()
print('after source filter:', lab.shape, '| pos', lab.is_positive.sum(), 'neg', (~lab.is_positive).sum())
print('positive sources:', lab[lab.is_positive].src.value_counts().to_dict())

# merge features
m=lab.merge(f,on=KEY,how='left',suffixes=('','_f')).merge(g,on=KEY,how='left')
gcols=['dnase_mean','rloop_mean','mappability_mean','atac_overlap']
fcols=[c for c in f.columns if c not in KEY+['region','center_trinuc'] and pd.api.types.is_numeric_dtype(f[c])]
print(f'\nf feature cols: {len(fcols)} | g feature cols: {len(gcols)}')
fcov=m[fcols].notna().all(axis=1).mean(); gcov=m[gcols].notna().all(axis=1).mean()
print(f'feature coverage: f={fcov:.2%} g={gcov:.2%}')

# keep rows with full features
ok=m[fcols].notna().all(axis=1) & m[gcols].notna().all(axis=1)
m=m[ok].copy()
print(f'rows with full f+g: {len(m)} (pos {m.is_positive.sum()} / neg {(~m.is_positive).sum()})')
print('by class:', m.groupby(['cls','is_positive']).size().to_dict())

m.to_parquet('data/processed/train_table_v1.parquet',index=False)
# also save the feature-column manifest
import json
json.dump({'f':fcols,'g':gcols},open('data/processed/feature_cols_v1.json','w'),indent=0)
print('\nwrote data/processed/train_table_v1.parquet + feature_cols_v1.json')
