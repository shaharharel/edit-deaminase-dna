"""Audit follow-up: does REPLICATION TIMING (Repli-seq, the mechanistically-correct magnitude covariate) predict
domain editing-burden RATE, where |RFD| did NOT? Closes the magnitude leg honestly. Poisson, uniform callable bg."""
import pandas as pd, numpy as np, pyBigWig
from scipy.stats import spearmanr
import statsmodels.api as sm
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','BE4_CT','BE4_tot','Parent_CT','Parent_tot'])
sp=sp[sp.motif=='TpC'].copy()
sp['callable']=(sp.BE4_tot>=10); sp['editor']=((sp.BE4_CT>=3)&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0))
sp['bin']=sp.pos//1_000_000
rt=pyBigWig.open('/tmp/poc_dna/hg19_hela_rt.bw'); CHR=set(rt.chroms())
dom=sp.groupby(['chrom','bin']).agg(n_editor=('editor','sum'),n_callable=('callable','sum'),
        cov=('BE4_tot','mean')).reset_index()
dom=dom[dom.n_callable>=50]
def binrt(ch,b):
    c2=ch if ch in CHR else ('chr'+ch if 'chr'+ch in CHR else None)
    if c2 is None: return np.nan
    s=int(b*1_000_000); e=min(s+1_000_000, rt.chroms(c2))
    try:
        v=rt.stats(c2,s,e,type="mean")[0]; return v
    except: return np.nan
dom['rt']=[binrt(c,b) for c,b in zip(dom.chrom,dom.bin)]
dom=dom.dropna(subset=['rt']); dom['rate']=dom.n_editor/dom.n_callable
print(f"=== DOMAN 1Mb domain-burden vs REPLICATION TIMING (Repli-seq, n_domains={len(dom)}) ===")
print("   (RT high = early-replicating; audit's correct magnitude covariate vs the null |RFD|)")
rho,pv=spearmanr(dom.rt,dom.rate)
print(f"  editor burden RATE vs RT: Spearman={rho:+.3f} p={pv:.2g}")
X=sm.add_constant(pd.DataFrame({'rt':(dom.rt-dom.rt.mean())/dom.rt.std(),
                                'logcov':np.log(dom['cov'].values)}))
m=sm.GLM(dom.n_editor,X,family=sm.families.Poisson(),offset=np.log(dom.n_callable)).fit()
print(f"  Poisson n_editor ~ RT(z)+logcov, offset log(callable): RT coef={m.params['rt']:+.3f} (p={m.pvalues['rt']:.2g}) IRR/SD={np.exp(m.params['rt']):.3f}")
print(f"    logcov coef={m.params['logcov']:+.3f} (p={m.pvalues['logcov']:.1e})")
dom['t']=pd.qcut(dom.rt,3,labels=['late','mid','early'])
strat=dom.groupby('t',observed=True).rate.mean()
print(f"  burden rate by RT tertile: {dict(strat.round(6))}  early/late fold={strat.get('early')/max(strat.get('late'),1e-12):.2f}x")
print("\n  READ: if RT coef sig>0 (early-rep -> more editing) after coverage control -> magnitude leg CLOSES (mechanism drives burden via timing);")
print("        if null -> mechanism is STRAND-ONLY (bounds claim honestly, consistent with GOTI RT-null).")
