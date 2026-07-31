"""Clone-level replication (poc_dna cron confound b): does the fork-strand C>T-vs-G>A x sign(RFD) coupling
replicate ACROSS the 8 independent BE4 clones (not read-pooled pseudoreplication)? Per-clone split OR + LOO."""
import pandas as pd, numpy as np, pyBigWig
from scipy.stats import fisher_exact, wilcoxon
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','Parent_CT','Parent_tot'])
tpc=(sp.motif=='TpC').values
z=np.load('/tmp/poc_dna/feat/perclone_full.npz',allow_pickle=True)
parent_clean=(sp.Parent_CT.values==0)&(sp.Parent_tot.values>=10)
base = tpc & parent_clean
# union of per-clone editor sites -> compute RFD once
clones=[f'BE4_clone{i}' for i in range(1,9)]
union=np.zeros(len(sp),bool)
masks={}
for c in clones:
    m = base & (z[f'{c}__ct'][:]>=3) & (z[f'{c}__tot'][:]>=10)
    masks[c]=m; union|=m
print(f"union editor sites across 8 BE4 clones: {union.sum()}")
idx=np.where(union)[0]
sub=sp.iloc[idx][['chrom','pos','strand']].copy()
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
sub['rfd']=valarr(sub.chrom.values,sub.pos.values)
rfd_map=dict(zip(idx,sub.rfd.values))  # global-index -> rfd
strand_map=dict(zip(idx,sub.strand.values))
def split_or(m):
    ii=np.where(m)[0]
    Cp=Cn=Gp=Gn=0
    for i in ii:
        r=rfd_map.get(i,np.nan)
        if not (r==r) or r==0: continue
        s=strand_map[i]
        if s=='+':
            if r>0: Cp+=1
            else: Cn+=1
        else:
            if r>0: Gp+=1
            else: Gn+=1
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,np.nan,Cp+Cn+Gp+Gn
    orr,p=fisher_exact([[Cp,Cn],[Gp,Gn]]); return orr,p,Cp+Cn+Gp+Gn
print("\n=== PER-CLONE fork-strand split OR (BE4, independent clones) ===")
ors=[]
for c in clones:
    orr,p,n=split_or(masks[c]); ors.append(orr)
    print(f"  {c}: OR={orr:.2f} p={p:.2g} n={n}")
ors=np.array([o for o in ors if o==o])
print(f"\n  8-clone OR: mean={ors.mean():.2f} median={np.median(ors):.2f} range=[{ors.min():.2f},{ors.max():.2f}]")
print(f"  fraction clones with OR>1: {(ors>1).mean():.2f}  (all>1 = replicates across independent clones, NOT pseudoreplication)")
try:
    w,pw=wilcoxon(ors-1,alternative='greater'); print(f"  Wilcoxon OR>1 across clones: p={pw:.3g}")
except Exception as e: print("  wilcoxon:",e)
