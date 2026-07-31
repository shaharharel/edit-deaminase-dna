import pandas as pd, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
d=pd.read_parquet('/tmp/poc_dna/goti/goti_rfd.parquet').dropna(subset=['rfd']).reset_index(drop=True)
tri=pd.get_dummies(d['tri'],prefix='tri').astype('float32')
CALL=['gc200','rep','in_gene','log_dist_tss','dnase']
y=d['label'].values.astype(int); chrom=d['Chrom'].values
def locho(cols):
    X=pd.concat([d[CALL+cols].astype('float32'),tri],axis=1).values
    oof=np.full(len(y),np.nan)
    for c in np.unique(chrom):
        tr=chrom!=c; te=chrom==c
        if y[tr].sum()==0: continue
        m=HistGradientBoostingClassifier(max_depth=3,learning_rate=0.05,max_iter=200,
              l2_regularization=1.0,min_samples_leaf=30,random_state=0)
        m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
    ok=~np.isnan(oof); return roc_auc_score(y[ok],oof[ok])
base=locho([])
print(f"=== GOTI mechanism decomposition (LOChr-out AUC, delta over callability base) ===")
print(f"  base (callability+access):          {base:.4f}")
for name,cols in [('+ on_exposed (DIRECTIONAL)',['on_exposed']),
                  ('+ absrfd (MAGNITUDE only) ',['absrfd']),
                  ('+ rfd (signed)            ',['rfd']),
                  ('+ on_exposed + absrfd     ',['on_exposed','absrfd']),
                  ('+ all {rfd,on_exp,absrfd} ',['rfd','on_exposed','absrfd'])]:
    a=locho(cols); print(f"  {name}: {a:.4f}  (delta {a-base:+.4f})")
print("\n  INTERP: if on_exposed-only ~= all-delta and >> absrfd-only -> DIRECTIONAL (lagging-strand) mechanism confirmed.")
print("          if absrfd-only carries it -> MAGNITUDE (near init/term zones), directional story WRONG.")
