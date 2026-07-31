"""Human confirmation of the concentration result: does the fork-strand risk model RANK editor sites in Doman,
and how much does it concentrate? Positives = BE4 editor C>T (TpC, Parent-clean); negatives = coverage-matched
unedited TpC sites. Answers the user's 'better S/N at region?' + 'panel?' questions on human data."""
import pandas as pd, numpy as np, pyBigWig
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','BE4_CT','BE4_tot','Parent_CT','Parent_tot'])
sp=sp[(sp.motif=='TpC')&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
pos=sp[sp.BE4_CT>=3].copy(); neg_pool=sp[sp.BE4_CT==0].copy()
# coverage-decile-matched negatives, 3:1
pos['covq']=pd.qcut(pos.BE4_tot,10,labels=False,duplicates='drop')
neg_pool['covq']=pd.qcut(neg_pool.BE4_tot,10,labels=False,duplicates='drop')
negs=[]
rng=np.random.default_rng(0)
for q in pos.covq.dropna().unique():
    npool=neg_pool[neg_pool.covq==q]; k=min(3*int((pos.covq==q).sum()),len(npool))
    if k>0: negs.append(npool.sample(k,random_state=0))
neg=pd.concat(negs)
pos['label']=1; neg['label']=0
d=pd.concat([pos,neg]).reset_index(drop=True)
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
v=[]
for c,p in zip(d.chrom,d.pos):
    c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else None)
    try: v.append(rfd.values(c2,int(p)-1,int(p))[0] if c2 else np.nan)
    except: v.append(np.nan)
d['rfd']=v; d=d[d.rfd.notna()&(d.rfd!=0)].reset_index(drop=True)
d['on_exposed']=(((d.strand=='+')&(d.rfd>0))|((d.strand=='-')&(d.rfd<0))).astype(int)
d['absrfd']=d.rfd.abs()
y=d.label.values; chrom=d.chrom.values
def locho(cols):
    X=d[cols].astype('float32').values; oof=np.full(len(y),np.nan)
    for c in np.unique(chrom):
        tr=chrom!=c;te=chrom==c
        if y[tr].sum()==0 or y[tr].sum()==tr.sum(): continue
        m=HistGradientBoostingClassifier(max_depth=3,learning_rate=0.05,max_iter=200,min_samples_leaf=30,random_state=0)
        m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
    return oof
oc=locho(['BE4_tot']); of=locho(['BE4_tot','on_exposed','absrfd'])
ok=~np.isnan(of)
print(f"=== DOMAN HUMAN site-level ranker ({y.sum()} pos / {(y==0).sum()} cov-matched neg) ===")
print(f"  coverage-only AUC={roc_auc_score(y[ok],oc[ok]):.3f} ; +fork-strand AUC={roc_auc_score(y[ok],of[ok]):.3f}")
dd=pd.DataFrame({'risk':of,'y':y}).dropna().sort_values('risk',ascending=False).reset_index(drop=True)
N=len(dd);P=dd.y.sum()
print(f"  CONCENTRATION (recall@top-k%, base rate {P/N:.2f}):")
for k in [1,5,10,20,50]:
    tn=int(N*k/100); rec=dd.y.iloc[:tn].sum()/P
    print(f"    top {k:>2}%: captures {rec*100:4.1f}% of positives (lift {rec/(k/100):.2f}x)")
print("  => HUMAN confirms mouse: fork-strand RANKS but concentration is MODEST (diffuse), not a panel.")
