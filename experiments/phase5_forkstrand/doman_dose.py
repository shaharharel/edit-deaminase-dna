"""|RFD| dose-response in Doman HUMAN: does the C>T-vs-G>A split OR increase with fork-polarization strength (|RFD|)?
Monotone increase = dose-response = strong causal evidence the coupling is replication-fork-driven."""
import pandas as pd, numpy as np, pyBigWig
from scipy.stats import fisher_exact
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','BE4_CT','BE4_tot','Parent_CT','Parent_tot'])
sp=sp[sp.motif=='TpC']
be=sp[(sp.BE4_CT>=3)&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
v=[]
for c,p in zip(be.chrom,be.pos):
    c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else None)
    try: v.append(rfd.values(c2,int(p)-1,int(p))[0] if c2 else np.nan)
    except: v.append(np.nan)
be['rfd']=v; be=be[be.rfd.notna()&(be.rfd!=0)]; be['absrfd']=be.rfd.abs()
be['oe']=((be.strand=='+')&(be.rfd>0))|((be.strand=='-')&(be.rfd<0))
def splitOR(df):
    ctp=df[df.strand=='+']; gap=df[df.strand=='-']
    Cp=int((ctp.rfd>0).sum());Cn=int((ctp.rfd<0).sum());Gp=int((gap.rfd>0).sum());Gn=int((gap.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,np.nan
    orr,p=fisher_exact([[Cp,Cn],[Gp,Gn]]); return orr,p
be['tert']=pd.qcut(be.absrfd,3,labels=['low|RFD|','mid','high|RFD|'])
print("=== DOMAN HUMAN |RFD| dose-response (fork polarization -> coupling strength) ===")
for t in ['low|RFD|','mid','high|RFD|']:
    d=be[be.tert==t]; orr,p=splitOR(d)
    print(f"  {t:9} (n={len(d)}, mean|RFD|={d.absrfd.mean():.2f}): on_exposed={d.oe.mean():.3f}  split OR={orr:.2f} p={p:.2g}")
# quintiles for finer dose curve
be['q5']=pd.qcut(be.absrfd,5,labels=False)
print("\n  on_exposed by |RFD| quintile:", [round(be[be.q5==q].oe.mean(),3) for q in range(5)])
from scipy.stats import spearmanr
rho,pv=spearmanr(be.absrfd,be.oe.astype(int))
print(f"  Spearman(|RFD|, on_exposed)={rho:+.3f} p={pv:.2g}  (monotone rise = dose-response)")
