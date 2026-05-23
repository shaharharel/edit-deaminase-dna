"""Add iPSC + HSPC DNase per-bin means to bins_v2 -> bins_v3 (for symmetric cell-type test + clinical demo)."""
import pandas as pd, numpy as np, pyBigWig, os
df=pd.read_parquet('/mnt/data/dna_features/bins_1mb_v2.parquet')
BIN=1_000_000
new={'dnase_ipsc':'/mnt/data/encode_tracks/ipsc/DNase.bigWig','dnase_hspc':'/mnt/data/encode_tracks/hspc/DNase.bigWig'}
bw={k:pyBigWig.open(v) for k,v in new.items() if os.path.exists(v)}
for k in bw: df[k]=np.nan
for i,r in df.iterrows():
    s=int(r.bin)*BIN; e=s+BIN
    for k,h in bw.items():
        try:
            if r.chrom in h.chroms():
                v=h.stats(r.chrom,s,min(e,h.chroms()[r.chrom]),type='mean')[0]
                df.at[i,k]=float(v) if v is not None else np.nan
        except: pass
df.to_parquet('/mnt/data/dna_features/bins_1mb_v3.parquet')
print('bins_v3:',df.shape)
print('coverage:',{k:round(df[k].notna().mean(),2) for k in bw})
