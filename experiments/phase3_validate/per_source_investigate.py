"""Why are some sources chromatin-predictable and others not? Per-source: dispersion (Gini),
raw correlation with accessibility, motif, assay type. Defines what 'predictable off-target risk' means."""
import numpy as np, pandas as pd
df=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','rloop','atac']).reset_index(drop=True)
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values
    if np.std(a)==0 or np.std(b)==0: return np.nan
    return float(np.corrcoef(a,b)[0,1])
def gini(x):
    x=np.sort(np.asarray(x,float)); n=len(x); c=np.cumsum(x)
    return float((n+1-2*np.sum(c)/c[-1])/n) if c[-1]>0 else 0.0
ASSAY={'ec_Doman':'orthogonal-R-loop+WGS','ec_McGrath':'iPSC clonal WGS de-novo','ec_Yu':'WGS off-target',
       'ec_Lei':'Detect-seq (genome-wide)','ec_Chen':'tdCBE'}
MOTIF={'ec_Doman':0.81,'ec_McGrath':0.90,'ec_Yu':0.17,'ec_Lei':0.48,'ec_Chen':0.01}
PRED={'ec_Doman':0.22,'ec_McGrath':0.006,'ec_Yu':0.64,'ec_Lei':0.40,'ec_Chen':0.23}  # from g_model_v2
print(f"{'source':10s} {'assay':24s} {'n':>6s} {'TpC':>5s} {'Gini':>5s} {'top10%':>6s} {'r_DNase':>8s} {'r_ATAC':>7s} {'r_Rloop':>8s} {'r_eligC':>8s} {'gPred':>6s}")
for s in ['ec_Doman','ec_McGrath','ec_Yu','ec_Lei','ec_Chen']:
    n=int(df[s].sum())
    if n<50: continue
    g=gini(df[s]); t10=df.nlargest(int(0.1*len(df)),s)[s].sum()/max(df[s].sum(),1)
    print(f"  {s[3:]:8s} {ASSAY[s]:24s} {n:6d} {MOTIF[s]:5.2f} {g:5.2f} {t10:6.2f} {sp(df[s],df.dnase):8.3f} {sp(df[s],df.atac):7.3f} {sp(df[s],df.rloop):8.3f} {sp(df[s],df.eligible_C):8.3f} {PRED[s]:6.3f}")
print("\n[read] high r_DNase/r_ATAC + high Gini (clustered) => chromatin-predictable.")
print("Compare to motif(TpC) and gPred: is predictability about accessibility-clustering, not motif cleanliness?")
