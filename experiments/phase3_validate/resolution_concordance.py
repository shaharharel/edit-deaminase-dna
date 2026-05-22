"""QA P1-1: resolution-vs-concordance curve. At what spatial resolution does guide-independent
CBE off-target burden become REPRODUCIBLE across independent assays? Sets the honest unit for the
safety gate. Mean pairwise Spearman across guide-indep CBE sources, chrom-block bootstrap CI."""
import numpy as np, pandas as pd
RNG=np.random.default_rng(0)
def spearman(a,b):
    a=pd.Series(a).rank().values; b=pd.Series(b).rank().values
    if a.std()==0 or b.std()==0: return np.nan
    return float(np.corrcoef(a,b)[0,1])

df=pd.read_parquet('data/processed/dna_offtarget_sites_v2.parquet')
df['src']=df.source_paper.str.replace(r'_iPSC_.*','_iPSC',regex=True)
df['cls']=df.apply(lambda r:'CBE' if (r.ref,r.alt) in [('C','T'),('G','A')] else 'x',axis=1)
CH=[f'chr{i}' for i in range(1,23)]+['chrX']
df=df[(df.cls=='CBE')&df.chrom.isin(CH)]
SRC=['Doman_BE4_pilot','McGrath_2019_iPSC','Yu_2020','Chen_2023_tdCBE']  # guide-indep CBE, Lei dropped (bug)
df=df[df.src.isin(SRC)]
print("CBE guide-indep sites by source:", df.src.value_counts().to_dict())

def curve(binsize, label):
    if binsize is None: return
    df['bin']=df.chrom+'_'+(df.pos//binsize).astype(str)
    counts={s:df[df.src==s].groupby('bin').size() for s in SRC}
    allb=sorted(set().union(*[set(c.index) for c in counts.values()]))
    binchrom={b:b.rsplit('_',1)[0] for b in allb}
    # pairwise mean spearman
    rhos=[]
    for i in range(len(SRC)):
        for j in range(i+1,len(SRC)):
            a=counts[SRC[i]].reindex(allb).fillna(0); b=counts[SRC[j]].reindex(allb).fillna(0)
            r=spearman(a.values,b.values)
            if not np.isnan(r): rhos.append(r)
    mean_rho=np.mean(rhos)
    # chrom-block bootstrap of the mean pairwise rho
    chroms=list(set(binchrom.values())); boot=[]
    for _ in range(300):
        pick=set(RNG.choice(chroms,len(chroms),replace=True))
        bb=[b for b in allb if binchrom[b] in pick]
        if len(bb)<20: continue
        rr=[]
        for i in range(len(SRC)):
            for j in range(i+1,len(SRC)):
                a=counts[SRC[i]].reindex(bb).fillna(0); b=counts[SRC[j]].reindex(bb).fillna(0)
                r=spearman(a.values,b.values)
                if not np.isnan(r): rr.append(r)
        if rr: boot.append(np.mean(rr))
    lo,hi=np.percentile(boot,[2.5,97.5]) if boot else (np.nan,np.nan)
    print(f"  {label:>8}: mean pairwise Spearman = {mean_rho:+.3f}  [95% CI {lo:+.3f}, {hi:+.3f}]  (n_bins={len(allb)})")
    return mean_rho

print("\n=== resolution-vs-concordance (guide-indep CBE, Lei excluded) ===")
for bs,lab in [(10_000,'10kb'),(100_000,'100kb'),(1_000_000,'1Mb'),(5_000_000,'5Mb'),(10_000_000,'10Mb')]:
    curve(bs,lab)
print("\n[interpretation] crossover resolution = finest scale where CI excludes 0 = honest unit")
print("for calibrated burden prediction. Finer than that = report intrinsic prior only.")
