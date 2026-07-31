#!/usr/bin/env python
"""GOTI feasibility — QA-CORRECTED (ml-code-reviewer). The claim is "replication-fork-strand BEATS a
callability+accessibility baseline", so the null/CI must be on the DELTA (auc_full - auc_call), NOT auc_full alone.
 - Paired bootstrap CI on observed delta (resample OOF pairs).
 - Permutation DELTA null: per perm, retrain BOTH models on same within-chrom-shuffled labels, collect af-ac.
 - Report pooled AND macro (per-chrom) AUC.
 - DROPPED the 1Mb domain-posfrac analysis: confounded (matched-neg spatial density tracks positive density,
   so posfrac ~= local matching ratio, not a true rate). Proper burden needs a uniform callable background (not available for GOTI matched negs).
"""
import pandas as pd, numpy as np, sys
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
rng=np.random.default_rng(0)
d=pd.read_parquet('/tmp/poc_dna/goti/goti_rfd.parquet').dropna(subset=['rfd']).reset_index(drop=True)
tri=pd.get_dummies(d['tri'],prefix='tri').astype('float32')
CALL=['gc200','rep','in_gene','log_dist_tss','dnase']; STRAND=['rfd','on_exposed','absrfd']
Xcall=pd.concat([d[CALL].astype('float32'),tri],axis=1).values
Xfull=pd.concat([d[CALL+STRAND].astype('float32'),tri],axis=1).values
y=d['label'].values.astype(int); chrom=d['Chrom'].values
NPERM=200; MI=200

def locho(X,yy):
    oof=np.full(len(yy),np.nan)
    for c in np.unique(chrom):
        tr=chrom!=c; te=chrom==c
        if yy[tr].sum()==0 or yy[tr].sum()==tr.sum(): continue
        m=HistGradientBoostingClassifier(max_depth=3,learning_rate=0.05,max_iter=MI,
              l2_regularization=1.0,min_samples_leaf=30,random_state=0)
        m.fit(X[tr],yy[tr]); oof[te]=m.predict_proba(X[te])[:,1]
    return oof
def pooled_auc(oof,yy):
    ok=~np.isnan(oof); return roc_auc_score(yy[ok],oof[ok])
def macro_auc(oof,yy):
    a=[]
    for c in np.unique(chrom):
        m=(chrom==c)&(~np.isnan(oof))
        if yy[m].sum()>0 and yy[m].sum()<m.sum(): a.append(roc_auc_score(yy[m],oof[m]))
    return np.mean(a)

oof_c=locho(Xcall,y); oof_f=locho(Xfull,y)
auc_c=pooled_auc(oof_c,y); auc_f=pooled_auc(oof_f,y); delta=auc_f-auc_c
print(f"=== GOTI LOChr-out (guide-indep CBE, mm10; 1609 pos / {int((y==0).sum())} detectability-matched neg) ===",flush=True)
print(f"  callability+access baseline: pooled AUC={auc_c:.4f}  macro={macro_auc(oof_c,y):.4f}",flush=True)
print(f"  + fork-strand:               pooled AUC={auc_f:.4f}  macro={macro_auc(oof_f,y):.4f}",flush=True)
print(f"  DELTA (strand incremental) = {delta:+.4f}",flush=True)

# paired bootstrap CI on delta (resample sites, recompute both AUCs on same indices)
ok=~np.isnan(oof_c)&~np.isnan(oof_f); yc=y[ok]; oc=oof_c[ok]; of=oof_f[ok]; idx=np.arange(len(yc))
bd=[]
for _ in range(2000):
    b=rng.choice(idx,len(idx),replace=True)
    if yc[b].sum()==0 or yc[b].sum()==len(b): continue
    bd.append(roc_auc_score(yc[b],of[b])-roc_auc_score(yc[b],oc[b]))
bd=np.array(bd); lo,hi=np.quantile(bd,[.025,.975])
print(f"  paired bootstrap delta 95% CI = [{lo:+.4f},{hi:+.4f}]  p(delta<=0)={(bd<=0).mean():.4f}",flush=True)

# permutation DELTA null: retrain BOTH on within-chrom-shuffled labels
print(f"  running {NPERM}-perm DELTA null (retrains both models)...",flush=True)
pd_=[]
for i in range(NPERM):
    yp=y.copy()
    for c in np.unique(chrom):
        ix=np.where(chrom==c)[0]; yp[ix]=rng.permutation(yp[ix])
    dc=locho(Xcall,yp); df_=locho(Xfull,yp)
    pd_.append(pooled_auc(df_,y)-pooled_auc(dc,y))
    if (i+1)%50==0: print(f"    perm {i+1}/{NPERM}",flush=True)
pd_=np.array(pd_); p=(1+ (pd_>=delta).sum())/(1+NPERM)
print(f"  perm DELTA null: mean={pd_.mean():+.4f} 95th={np.quantile(pd_,0.95):+.4f}  p(perm_delta>=obs)={p:.4f}",flush=True)
verdict = "STRAND BEATS callability baseline (delta CI>0 AND perm p<0.05)" if (lo>0 and p<0.05) else \
          ("strand adds (CI>0) but perm-null marginal" if lo>0 else "strand does NOT beat callability")
print(f"\n  VERDICT: {verdict}",flush=True)
# calibration (honest caveat: to the 25% matched base rate, NOT absolute genomic prior)
dd=pd.DataFrame({'risk':oof_f,'label':y}).dropna()
dd['q']=pd.qcut(dd.risk,10,labels=False,duplicates='drop')
cal=dd.groupby('q').agg(pred=('risk','mean'),obs=('label','mean'),n=('label','size'))
ece=(cal.n*(cal.pred-cal.obs).abs()).sum()/cal.n.sum()
print(f"  calibration ECE={ece:.4f} (to 25% MATCHED base rate — NOT absolute genomic off-target prior)",flush=True)
print("DONE goti_domain_v2",flush=True)
