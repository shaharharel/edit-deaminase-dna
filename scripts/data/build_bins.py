"""g v2: per-1Mb bins + TpC motif count (control) + histone(.bigWig) + ATAC(bed) + repliseq(if present)."""
import pandas as pd, numpy as np, pyBigWig, pyfaidx
from collections import defaultdict
fa=pyfaidx.Fasta('/mnt/data/ref/hg38/hg38.fa')
BIN=1_000_000; CH=[f'chr{i}' for i in range(1,23)]+['chrX']
df=pd.read_parquet('/mnt/data/dna_features/canonical_sites.parquet')
GI=['Doman_BE4_pilot','McGrath_2019_iPSC','Yu_2020','Lei_2021','Chen_2023_tdCBE']
ce=df[(df.edit_class=='CBE')&(df.src.isin(GI))]
ec=defaultdict(int); ecs=defaultdict(lambda:defaultdict(int))
for r in ce.itertuples(): ec[(r.chrom,r.pos//BIN)]+=1; ecs[(r.chrom,r.pos//BIN)][r.src]+=1
import os
tracks={'rloop':'rloop/HEK293_MapR_comb.bw','dnase':'dnase/ENCFF148BGE.bigWig',
  'h3k27ac':'histone/H3K27ac.bigWig','h3k4me3':'histone/H3K4me3.bigWig',
  'phylop':'conservation/phyloP100way.bw','phastcons':'conservation/phastCons100way.bw',
  'mappability':'mappability/k100.Umap.MultiTrackMappability.bw'}
# repliseq if present
for f in os.listdir('/mnt/data/encode_tracks/repliseq') if os.path.isdir('/mnt/data/encode_tracks/repliseq') else []:
    if f.endswith(('.bw','.bigWig')): tracks['repliseq']='repliseq/'+f; break
bw={}
for k,v in tracks.items():
    p='/mnt/data/encode_tracks/'+v
    if os.path.exists(p):
        try: bw[k]=pyBigWig.open(p)
        except Exception as e: print('skip',k,e)
# ATAC bed -> per-bin covered fraction
atac=defaultdict(int)
abed='/mnt/data/encode_tracks/atac/ATC.Kid.20.AllAg.293.bed'
if os.path.exists(abed):
    for line in open(abed):
        p=line.split('\t')
        if len(p)<3 or not p[0].startswith('chr'): continue
        try: c,s,e=p[0],int(p[1]),int(p[2])
        except: continue
        atac[(c,s//BIN)]+=min(e,(s//BIN+1)*BIN)-s
rows=[]
for c in CH:
    L=len(fa[c])
    for b in range(L//BIN):
        s=b*BIN; e=min(s+BIN,L); seq=fa[c][s:e].seq.upper()
        nC=seq.count('C')+seq.count('G')
        if nC<10000: continue
        tpc=seq.count('TC')+seq.count('GA')   # TpC motif opportunity (both strands)
        feat={'chrom':c,'bin':b,'eligible_C':nC,'tpc_count':tpc,'edit_count':ec.get((c,b),0),
              'atac':atac.get((c,b),0)/BIN}
        for k,h in bw.items():
            try:
                v=h.stats(c,s,e,type='mean')[0] if c in h.chroms() else None
                feat[k]=float(v) if v is not None else np.nan
            except: feat[k]=np.nan
        for srx in GI: feat['ec_'+srx.split('_')[0]]=ecs[(c,b)].get(srx,0)
        rows.append(feat)
out=pd.DataFrame(rows); out.to_parquet('/mnt/data/dna_features/bins_1mb_v2.parquet')
print("bins:",out.shape,"edits:",out.edit_count.sum())
print("features:",[c for c in out.columns if c not in ('chrom','bin')][:20])
print("track cov:", {k:round(out[k].notna().mean(),2) for k in bw})
print("atac>0 bins:", int((out.atac>0).sum()))
