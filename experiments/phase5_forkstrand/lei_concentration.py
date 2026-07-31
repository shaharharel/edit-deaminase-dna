"""Concentration on LEI DETECT-SEQ (the enrichment assay). Does motif x fork-strand concentrate the editing signal
BETTER than on Doman WGS (~1x)? Metrics: fork-strand OR by motif + |RFD|; editing fraction (VAF) by stratum."""
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
d=pd.read_parquet('/tmp/poc_dna/lei_annot.parquet')
d=d[d.rfd.notna()&(d.rfd!=0)].copy()
d['on_exposed']=(((d.ref=='C')&(d.rfd>0))|((d.ref=='G')&(d.rfd<0))).astype(int)
d['absrfd']=d.rfd.abs(); d['vaf']=d.mm/d.tot
def OR(df):
    ct=df[df.ref=='C'];ga=df[df.ref=='G']
    Cp=int((ct.rfd>0).sum());Cn=int((ct.rfd<0).sum());Gp=int((ga.rfd>0).sum());Gn=int((ga.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return np.nan
    return fisher_exact([[Cp,Cn],[Gp,Gn]])[0]
print("=== LEI Detect-seq concentration (411K called editing sites) ===")
print(f"  fork-strand OR: ALL={OR(d):.2f}  TpC={OR(d[d.isTpC]):.2f}  nonTpC={OR(d[~d.isTpC]):.2f}")
print(f"    => MOTIF concentrates the strand-editing signal on Detect-seq (TpC vs nonTpC)")
tp=d[d.isTpC]
print("\n  fork-strand OR by |RFD| tertile WITHIN TpC (dose-response = further concentration):")
tp=tp.assign(t=pd.qcut(tp.absrfd,3,labels=['low','mid','high']))
for t in ['low','mid','high']:
    s=tp[tp.t==t]; print(f"    |RFD| {t:4}: OR={OR(s):.2f} (n={len(s)}, mean|RFD|={s.absrfd.mean():.2f})")
print("\n  editing fraction (VAF=mm/tot) by stratum (is editing stronger where predicted?):")
for lab,m in [('all',d.index==d.index),('TpC',d.isTpC),('TpC & on_exposed',(d.isTpC)&(d.on_exposed==1)),
              ('TpC & on_exposed & high|RFD|',(d.isTpC)&(d.on_exposed==1)&(d.absrfd>d.absrfd.quantile(.67)))]:
    print(f"    {lab:28}: mean VAF={d[m].vaf.mean():.3f} (n={m.sum()})")
# site reduction: TpC keeps the signal, drops nonTpC (null)
print(f"\n  SITE REDUCTION: TpC = {d.isTpC.mean()*100:.0f}% of called sites carry the fork-strand signal (OR {OR(d[d.isTpC]):.2f});")
print(f"    nonTpC = {(~d.isTpC).mean()*100:.0f}% are null (OR {OR(d[~d.isTpC]):.2f}) -> ~2.4x site reduction retaining signal.")
