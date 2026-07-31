"""DECISIVE QA control: is the Lei OR=1.25 real editing or germline-strand-asymmetry confound?
Germline het SNPs sit at VAF~0.5 (homs ~1.0); guide-independent editor edits are sub-clonal LOW-VAF.
Re-run C>T-vs-G>A x sign(RFD) split within VAF bins. If OR lives at VAF~0.5 and collapses at low-VAF -> GERMLINE CONFOUND."""
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
d=pd.read_parquet('/tmp/poc_dna/lei/lei_editor_sites.parquet')
d=d[d.rfd.notna()&(d.rfd!=0)].copy()
d['vaf']=d.mm/d.tot
print(f"=== VAF distribution of {len(d)} 'positives' ===")
print("  VAF histogram (bins):")
h,edges=np.histogram(d.vaf,bins=[0,0.05,0.1,0.15,0.2,0.35,0.5,0.65,0.8,1.01])
for i in range(len(h)):
    print(f"    [{edges[i]:.2f},{edges[i+1]:.2f}): {h[i]:>8,} ({100*h[i]/len(d):.1f}%)")
print(f"  median VAF={d.vaf.median():.3f} ; frac VAF>0.35 (germline-like)={ (d.vaf>0.35).mean():.3f} ; frac VAF<0.15 (editor-like)={(d.vaf<0.15).mean():.3f}")
def splitOR(df):
    ct=df[df.ref=='C']; ga=df[df.ref=='G']
    Cp=int((ct.rfd>0).sum());Cn=int((ct.rfd<0).sum());Gp=int((ga.rfd>0).sum());Gn=int((ga.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,np.nan,len(df)
    orr,p=fisher_exact([[Cp,Cn],[Gp,Gn]]); return orr,p,len(df)
print("\n=== C>T-vs-G>A x sign(RFD) split BY VAF BIN ===")
bins=[(0,0.10,'editor-like <0.10'),(0.10,0.20,'low 0.10-0.20'),(0.20,0.35,'mid 0.20-0.35'),
      (0.35,0.65,'GERMLINE-het ~0.5'),(0.65,1.01,'GERMLINE-hom ~1.0')]
for lo,hi,lab in bins:
    sub=d[(d.vaf>=lo)&(d.vaf<hi)]; orr,p,n=splitOR(sub)
    print(f"  VAF [{lo:.2f},{hi:.2f}) {lab:20}: n={n:>8,} OR={orr:.3f} p={p:.2g}")
print("\nREAD: if OR>1 concentrated at VAF~0.5/1.0 and ~1.0 at low-VAF -> GERMLINE confound (downgrade Lei).")
print("      if OR>1 at low-VAF (editor-like) -> real editing signal survives.")
