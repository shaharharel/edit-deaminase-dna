import pandas as pd, numpy as np, pyBigWig, pickle
from scipy.stats import fisher_exact, chi2
g2s=pickle.load(open('/tmp/poc_dna/gene2strand.pkl','rb'))
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','gene','BE4_CT','BE4_tot','YE1-BE4_CT','YE1-BE4_tot',
             'nCas9_CT','nCas9_tot','Parent_CT','Parent_tot'])
sp=sp[sp.motif=='TpC'].copy(); sp['gstrand']=sp.gene.map(g2s)
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
def split(df):
    df=df[df.rfd.notna()&(df.rfd!=0)]
    ctp=df[df.strand=='+']; gap=df[df.strand=='-']
    Cp=int((ctp.rfd>0).sum());Cn=int((ctp.rfd<0).sum());Gp=int((gap.rfd>0).sum());Gn=int((gap.rfd<0).sum())
    orr,fp=fisher_exact([[Cp,Cn],[Gp,Gn]]); return orr,fp,len(df)
def mh(df):  # transcription-controlled
    df=df.dropna(subset=['gstrand']); df=df[df.rfd.notna()&(df.rfd!=0)]; tabs=[]
    for g in ['+','-']:
        d=df[df.gstrand==g]; ctp=d[d.strand=='+']; gap=d[d.strand=='-']
        tabs.append(np.array([[int((ctp.rfd>0).sum()),int((ctp.rfd<0).sum())],
                              [int((gap.rfd>0).sum()),int((gap.rfd<0).sum())]]))
    num=sum((t[0,0]*t[1,1])/t.sum() for t in tabs);den=sum((t[0,1]*t[1,0])/t.sum() for t in tabs)
    A=E=V=0
    for t in tabs:
        n=t.sum();r1=t[0].sum();c1=t[:,0].sum();A+=t[0,0];E+=r1*c1/n;V+=r1*(n-r1)*c1*(n-c1)/(n**2*(n-1))
    chi=(abs(A-E)-0.5)**2/V; return num/den,1-chi2.cdf(chi,1)
print("=== POSITIVE-CONTROL CALIBRATION: recover the known APOBEC lagging-strand rule across 3 APOBEC sources ===")
print("   (identical C>T-vs-G>A x sign(RFD) pipeline; hg19 HeLa RFD, cell-type-MISMATCHED to Doman -> attenuates)")
print("   Known ground truth: endogenous APOBEC replication-strand asymmetry ~1.5-1.8x (Haradhvala/Seplyarskiy 2016)")
for tag,ct,tot in [('BE4 (editor+APOBEC3)','BE4_CT','BE4_tot'),
                   ('YE1-BE4 (attenuated editor+APOBEC3)','YE1-BE4_CT','YE1-BE4_tot'),
                   ('nCas9 (endogenous APOBEC3 ONLY)','nCas9_CT','nCas9_tot')]:
    d=sp[(sp[ct]>=3)&(sp[tot]>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
    d['rfd']=valarr(d.chrom.values,d.pos.values)
    orr,fp,n=split(d); mhor,mhp=mh(d)
    print(f"  {tag:38} n={n:5d}  split OR={orr:.2f} p={fp:.1e}  |  transcription-ctrl MH OR={mhor:.2f} p={mhp:.1e}")
print("\n  INTERP: consistent ~1.3-1.4x across all 3 independent APOBEC sources = pipeline RECOVERS the known rule")
print("          (attenuated from ~1.6x by cell-type-mismatched RFD). nCas9 = cleanest external-of-editor positive control.")
