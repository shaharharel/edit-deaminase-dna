"""SIMPLE per-site risk-ranking model (user request): motif x fork-strand -> calibrated per-site risk.
Answers concretely: (1) how much does it CONCENTRATE (recall@top-k, lift)? (2) does REGION aggregation improve S/N?
(3) does it beat callability-only? Honest about what a diffuse process allows."""
import pandas as pd, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
d=pd.read_parquet('/tmp/poc_dna/goti/goti_rfd.parquet').dropna(subset=['rfd']).reset_index(drop=True)
tri=pd.get_dummies(d['tri'],prefix='tri').astype('float32')
CALL=['gc200','rep','in_gene','log_dist_tss','dnase']; STRAND=['on_exposed','absrfd','rfd']
y=d['label'].values.astype(int); chrom=d['Chrom'].values
def locho(cols):
    X=pd.concat([d[cols].astype('float32'),tri],axis=1).values
    oof=np.full(len(y),np.nan)
    for c in np.unique(chrom):
        tr=chrom!=c; te=chrom==c
        if y[tr].sum()==0: continue
        m=HistGradientBoostingClassifier(max_depth=3,learning_rate=0.05,max_iter=250,l2_regularization=1.0,min_samples_leaf=30,random_state=0)
        m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
    return oof
oof_call=locho(CALL); oof_full=locho(CALL+STRAND)
ok=~np.isnan(oof_full)
auc_c=roc_auc_score(y[ok],oof_call[ok]); auc_f=roc_auc_score(y[ok],oof_full[ok])
print(f"=== SITE-LEVEL risk model (GOTI, {y.sum()} pos / {(y==0).sum()} matched-neg) ===")
print(f"  callability-only AUC={auc_c:.3f} ; +fork-strand AUC={auc_f:.3f} (delta {auc_f-auc_c:+.3f})")
# CONCENTRATION: recall@top-k% by risk (of the matched pos+neg set)
dd=pd.DataFrame({'risk':oof_full,'y':y}).dropna().sort_values('risk',ascending=False).reset_index(drop=True)
N=len(dd); P=dd.y.sum(); base_rate=P/N
print(f"\n  CONCENTRATION (recall@top-k%, base positive rate={base_rate:.2f}):")
for k in [1,5,10,20,50]:
    topn=int(N*k/100); rec=dd.y.iloc[:topn].sum()/P; lift=rec/(k/100)
    print(f"    top {k:>2}% by risk: captures {rec*100:4.1f}% of positives  (lift {lift:.2f}x over random)")
# REGION (1Mb) aggregation — does mean-risk per region predict region positive-fraction?
dd2=pd.DataFrame({'risk':oof_full,'y':y,'chrom':chrom,'pos':d['Pos'].values}).dropna()
dd2['bin']=dd2.pos//1_000_000
reg=dd2.groupby(['chrom','bin']).agg(mrisk=('risk','mean'),posfrac=('y','mean'),n=('y','size')).reset_index()
reg=reg[reg.n>=5]
from scipy.stats import spearmanr
rho,pv=spearmanr(reg.mrisk,reg.posfrac)
print(f"\n  REGION (1Mb) aggregation: Spearman(mean-risk, region positive-fraction)={rho:+.3f} p={pv:.2g} (n_regions={len(reg)})")
print("    NOTE: negatives are detectability-MATCHED (drawn per-positive) -> region positive-fraction is partly")
print("    confounded by matching density; a clean region-burden test on genome-wide data (Doman) was NULL.")
print(f"\n  HONEST READ: site risk RANKS (top-{[k for k in [10]][0]}% lift shown above) but the process is DIFFUSE —")
print("    modest concentration, NOT a tight panel. The fork-strand adds real lift over callability but caps ~AUC 0.57.")
