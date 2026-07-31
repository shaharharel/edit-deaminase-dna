import pandas as pd, numpy as np, pyBigWig, pickle
from scipy.stats import fisher_exact, chi2
g2s=pickle.load(open('/tmp/poc_dna/gene2strand.pkl','rb'))
print("OR4F5 strand (should be -):", g2s.get('OR4F5'), "| SAMD11 (should be +):", g2s.get('SAMD11'))
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','gene','motif','strand','BE4_CT','BE4_tot','nCas9_CT','Parent_CT','Parent_tot'])
sp=sp[sp.motif=='TpC']; sp['gstrand']=sp.gene.map(g2s)
cov=sp.gstrand.notna().mean(); print(f"gene-strand coverage of sites: {cov:.3f}")
be=sp[(sp.BE4_CT>=3)&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].dropna(subset=['gstrand']).copy()
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
def valarr(chrom,pos):
    out=np.full(len(pos),np.nan)
    for i,(c,p) in enumerate(zip(chrom,pos)):
        c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else None)
        if c2 is None: continue
        try:
            v=rfd.values(c2,int(p)-1,int(p))[0]
            if v==v: out[i]=v
        except: pass
    return out
be['rfd']=valarr(be.chrom.values,be.pos.values); be=be[be.rfd.notna()&(be.rfd!=0)]
print(f"\nBE4 editor sites n={len(be)}  +genes:{(be.gstrand=='+').sum()} -genes:{(be.gstrand=='-').sum()}")
# sanity: is genomic C-strand correlated with gene strand? (transcription confound magnitude)
print("corr(genomic C-strand==+ , gene-strand==+):",
      round(np.corrcoef((be.strand=='+').astype(int),(be.gstrand=='+').astype(int))[0,1],3))

def split_or(df):
    ctp=df[df.strand=='+']; gap=df[df.strand=='-']
    Cp=int((ctp.rfd>0).sum());Cn=int((ctp.rfd<0).sum());Gp=int((gap.rfd>0).sum());Gn=int((gap.rfd<0).sum())
    t=np.array([[Cp,Cn],[Gp,Gn]])
    if min(Cp+Cn,Gp+Gn)==0: return np.nan,np.nan,t
    orr,fp=fisher_exact(t); return orr,fp,t
print("\n=== MH: C>T-vs-G>A x sign(RFD) STRATIFIED by TRUE gene strand (transcription control) ===")
tabs=[]
for gstr in ['+','-']:
    d=be[be.gstrand==gstr]; orr,fp,t=split_or(d)
    print(f"  {gstr}-gene: n={len(d)} 2x2{t.tolist()} OR={orr:.2f} p={fp:.3g}")
    tabs.append(t)
num=sum((t[0,0]*t[1,1])/t.sum() for t in tabs); den=sum((t[0,1]*t[1,0])/t.sum() for t in tabs)
mh_or=num/den
A=E=V=0
for t in tabs:
    n=t.sum();r1=t[0].sum();c1=t[:,0].sum();A+=t[0,0];E+=r1*c1/n;V+=r1*(n-r1)*c1*(n-c1)/(n**2*(n-1))
chi=(abs(A-E)-0.5)**2/V; p_cmh=1-chi2.cdf(chi,1)
print(f"\n  MH pooled OR={mh_or:.2f}  CMH chi2={chi:.2f} p={p_cmh:.3g}")
print("  VERDICT:", "REPLICATION-FORK SURVIVES transcription control -> lagging-strand mechanism supported"
      if (mh_or>1.15 and p_cmh<0.05) else "collapses -> transcription-coupled; keep weaker 'RFD predicts location'")
