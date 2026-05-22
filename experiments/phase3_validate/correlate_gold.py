"""P3 headline: does predicted INTRINSIC (sequence-only) per-gene CBE susceptibility correlate
with the REAL measured BE4 WGS editing index? With QA-mandated controls:
  - chromosome-block bootstrap CI (spatial autocorrelation)
  - permutation null (shuffle gene labels)
  - report at honest resolution; filter by gold coverage."""
import numpy as np, pandas as pd
RNG=np.random.default_rng(0)

pred=pd.read_parquet('data/processed/pred_gene_burden_v1.parquet')
gold=pd.read_parquet('data/processed/be4_gold_index.parquet')
# gene -> chrom map (for block bootstrap)
rg=pd.read_csv('/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt',sep='\t',header=None,
   names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','name2','s1','s2','frames'])
g2c=rg.groupby('name2').chrom.first()

def spearman(a,b):
    a=pd.Series(a).rank().values; b=pd.Series(b).rank().values
    return float(np.corrcoef(a,b)[0,1])

m=pred.merge(gold,on='gene',how='inner')
m['chrom']=m.gene.map(g2c)
m=m[m.chrom.notna()]
print(f"genes overlapping pred+gold: {len(m)}")
# coverage filter: enough WGS coverage to trust the index
for mincov in [100, 1000, 5000]:
    s=m[m.cov_treated>=mincov]
    if len(s)<50: continue
    rho=spearman(s.pred_burden, s.excess_index_pct)
    rho_t=spearman(s.pred_burden, s.index_treated_pct)
    # permutation null
    null=[spearman(s.pred_burden.values, RNG.permutation(s.excess_index_pct.values)) for _ in range(1000)]
    p=(np.sum(np.abs(null)>=abs(rho))+1)/1001
    # chromosome-block bootstrap CI
    chroms=s.chrom.unique(); boot=[]
    for _ in range(1000):
        pick=RNG.choice(chroms,len(chroms),replace=True)
        bs=pd.concat([s[s.chrom==c] for c in pick])
        if bs.pred_burden.nunique()>2: boot.append(spearman(bs.pred_burden,bs.excess_index_pct))
    lo,hi=np.percentile(boot,[2.5,97.5])
    print(f"\n[cov>={mincov}] n={len(s)}  Spearman(pred, EXCESS index)={rho:+.3f}  [95% CI {lo:+.3f},{hi:+.3f}]  perm-null p={p:.3f}")
    print(f"            Spearman(pred, treated index)={rho_t:+.3f}")
    # GC confound control (QA P0-1): does pred add beyond GC?
    if 'gc_mean' in s.columns:
        rho_gc=spearman(s.gc_mean, s.excess_index_pct)
        def _resid(y,x):
            x=pd.Series(x).rank().values; y=pd.Series(y).rank().values
            b=np.polyfit(x,y,1); return y-(b[0]*x+b[1])
        rp=_resid(s.pred_burden.values,s.gc_mean.values); rgd=_resid(s.excess_index_pct.values,s.gc_mean.values)
        print(f"            GC-only vs gold={rho_gc:+.3f} | pred-vs-gold PARTIAL(|GC)={float(np.corrcoef(rp,rgd)[0,1]):+.3f}")
    # top-decile enrichment: are high-pred genes enriched in high-gold?
    s=s.copy(); s['pred_top']=s.pred_burden>=s.pred_burden.quantile(0.9)
    hi_gold=s.excess_index_pct>=s.excess_index_pct.quantile(0.9)
    from math import comb
    ov=int((s.pred_top&hi_gold).sum()); K=int(s.pred_top.sum()); N=len(s); n=int(hi_gold.sum())
    exp=K*n/N
    print(f"            top-10% pred vs top-10% gold overlap={ov} (expected {exp:.1f}, {ov/max(exp,1):.1f}x)")
print("\n[interpretation] positive rho with CI excluding 0 + small perm-p => sequence-intrinsic")
print("susceptibility predicts real BE4 burden. Modest rho expected (f only; g not yet added).")
