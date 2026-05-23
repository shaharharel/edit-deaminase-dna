"""Rigorous g (per scientific-analyst): Poisson GLM + eligible-TpC offset, pool Yu+Lei (bulk only),
within-source noise ceiling, delta-over-baseline + block-bootstrap CI, circular-shift null, leave-one-source-out."""
import numpy as np, pandas as pd
from bisect import bisect_right
RNG=np.random.default_rng(0)
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values
    return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['obs']=b.ec_Yu+b.ec_Lei            # pooled BULK only (drop McGrath clonal, Doman orthogonal a-priori)
b['off']=np.log(b.tpc_count+1)       # opportunity offset
def poisson_glm(X,y,off,it=600,lr=0.1,l2=1e-2):
    X=(X-X.mean(0))/(X.std(0)+1e-6); n,d=X.shape; w=np.zeros(d); a=np.log(y.mean()+1e-6)
    for _ in range(it):
        mu=np.exp(a+off+X@w); g=mu-y
        w-=lr*(X.T@g/n+l2*w); a-=lr*g.mean()
    return a,w,X.mean(0) if False else (a,w)
def loco_pred(feats):
    pred=np.zeros(len(b))
    for ho in b.chrom.unique():
        tr=b[b.chrom!=ho];te=b[b.chrom==ho]
        mu_tr=tr[feats].mean(0).values; sd=tr[feats].std(0).values+1e-6
        Xtr=((tr[feats].values-mu_tr)/sd); Xte=((te[feats].values-mu_tr)/sd)
        a=np.log(tr.obs.mean()+1e-6); w=np.zeros(len(feats))
        for _ in range(600):
            mu=np.exp(a+tr.off.values+Xtr@w); g=mu-tr.obs.values
            w-=0.1*(Xtr.T@g/len(tr)+1e-2*w); a-=0.1*g.mean()
        pred[te.index]=a+te.off.values+Xte@w
    return pred
CHROM=['dnase','atac','rloop','h3k27ac']; BASE=['mappability']
# baselines
p_off=b.off.values
p_base=loco_pred(BASE)            # offset+mappability
p_chr=loco_pred(BASE+CHROM)       # +chromatin
r_off=sp(p_off,b.obs); r_base=sp(p_base,b.obs); r_chr=sp(p_chr,b.obs)
print(f"\n=== Poisson GLM LOCO Spearman (target=Yu+Lei pooled, offset=TpC) ===")
print(f"  opportunity offset only : {r_off:.3f}")
print(f"  + mappability           : {r_base:.3f}")
print(f"  + chromatin             : {r_chr:.3f}   (delta over base = {r_chr-r_base:+.3f})")
# block bootstrap CI on chromatin model + delta
boot=[];bootd=[]
chs=b.chrom.unique()
for _ in range(500):
    pick=RNG.choice(chs,len(chs),replace=True)
    idx=np.concatenate([np.where(b.chrom.values==c)[0] for c in pick])
    boot.append(sp(p_chr[idx],b.obs.values[idx])); bootd.append(sp(p_chr[idx],b.obs.values[idx])-sp(p_base[idx],b.obs.values[idx]))
print(f"  chromatin 95% CI: [{np.percentile(boot,2.5):.3f},{np.percentile(boot,97.5):.3f}]  delta 95% CI: [{np.percentile(bootd,2.5):+.3f},{np.percentile(bootd,97.5):+.3f}]")
# circular-shift null (spatial autocorrelation control)
null=[]
for _ in range(200):
    sh=int(RNG.integers(50,len(b)-50)); shifted=np.roll(b.obs.values,sh)
    null.append(sp(p_chr,shifted))
print(f"  circular-shift NULL: mean={np.mean(null):+.3f} 95th={np.percentile(null,95):+.3f}  (real {r_chr:.3f} should be far above)")
# within-source split-half NOISE CEILING (Yu)
sites=pd.read_parquet('data/processed/canonical_sites.parquet')
yu=sites[(sites.src=='Yu_2020')&(sites.edit_class=='CBE')].copy(); yu['k']=yu.chrom+'_'+(yu.pos//1_000_000).astype(str)
h=RNG.random(len(yu))<0.5
ca=yu[h].groupby('k').size(); cb=yu[~h].groupby('k').size()
allk=sorted(set(ca.index)|set(cb.index))
print(f"\n  NOISE CEILING (Yu split-half bin reliability): {sp(ca.reindex(allk).fillna(0),cb.reindex(allk).fillna(0)):.3f}")
# leave-one-source-out: fit on Yu bins, predict Lei
def fit_one(tgt):
    mu=b[CHROM+BASE].mean(0).values;sd=b[CHROM+BASE].std(0).values+1e-6;X=(b[CHROM+BASE].values-mu)/sd
    a=np.log(b[tgt].mean()+1e-6);w=np.zeros(len(CHROM+BASE));off=np.log(b.tpc_count.values+1)
    for _ in range(600):
        m=np.exp(a+off+X@w);g=m-b[tgt].values;w-=0.1*(X.T@g/len(b)+1e-2*w);a-=0.1*g.mean()
    return a+off+X@w
print(f"  leave-one-source-out: fit Yu -> predict Lei: {sp(fit_one('ec_Yu'),b.ec_Lei):.3f} | fit Lei -> predict Yu: {sp(fit_one('ec_Lei'),b.ec_Yu):.3f}")
