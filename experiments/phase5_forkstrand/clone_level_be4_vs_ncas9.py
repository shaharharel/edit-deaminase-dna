"""DECISIVE deaminase-specificity control (QA scientific-analyst endorsed).
Does the fork-strand C>T-vs-G>A x sign(RFD) coupling differ between BE4 (deaminase + nCas9) and nCas9-only
(nCas9, NO deaminase) clones? This isolates the DEAMINASE (Cas9 held constant) — unlike Lei NT-vs-Empty which
confounds deaminase with Cas9/vector/culture. Clone-level (8 BE4 vs 7 nCas9), per-clone OR then unweighted
Mann-Whitney (NOT read-pooled pseudoreplication, NOT variance-weighted).
Ascertainment matched: same coverage filter (ct>=3,tot>=10), Parent-clean, TpC, same RFD. Composition null in this
genic universe is ~1.02 flat (identical for both cohorts) so raw per-clone ORs are already composition-clean."""
import pandas as pd, numpy as np, pyBigWig
from scipy.stats import fisher_exact, mannwhitneyu
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','Parent_CT','Parent_tot'])
tpc=(sp.motif=='TpC').values
z=np.load('/tmp/poc_dna/feat/perclone_ct.npz',allow_pickle=True)
parent_clean=(sp.Parent_CT.values==0)&(sp.Parent_tot.values>=10)
base=tpc&parent_clean
BE4=[f'BE4_clone{i}' for i in range(1,9)]
NC =[f'nCas9_clone{i}' for i in range(1,8)]
masks={}; union=np.zeros(len(sp),bool)
for c in BE4+NC:
    m=base&(z[f'{c}__ct'][:]>=3)&(z[f'{c}__tot'][:]>=10)
    masks[c]=m; union|=m
idx=np.where(union)[0]
print(f"union editor sites (8 BE4 + 7 nCas9): {union.sum():,}")
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
rfd_map=dict(zip(idx,sub.rfd.values)); strand_map=dict(zip(idx,sub.strand.values))
def split_or(m):
    Cp=Cn=Gp=Gn=0
    for i in np.where(m)[0]:
        r=rfd_map.get(i,np.nan)
        if not (r==r) or r==0: continue
        s=strand_map[i]
        if s=='+':
            if r>0: Cp+=1
            else: Cn+=1
        else:
            if r>0: Gp+=1
            else: Gn+=1
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,Cp+Cn+Gp+Gn
    return fisher_exact([[Cp,Cn],[Gp,Gn]])[0],Cp+Cn+Gp+Gn
def cohort(names,tag):
    print(f"\n=== {tag} per-clone strand OR ===")
    ors=[]
    for c in names:
        o,n=split_or(masks[c])
        print(f"  {c}: OR={o:.3f} n={n}");
        if o==o: ors.append(o)
    return np.array(ors)
be=cohort(BE4,'BE4 (deaminase+nCas9)'); nc=cohort(NC,'nCas9-only (NO deaminase)')
print(f"\n=== DEAMINASE-SPECIFICITY (clone-level, unweighted) ===")
print(f"  BE4   n={len(be)} median OR={np.median(be):.3f} mean={be.mean():.3f} [{be.min():.2f},{be.max():.2f}]")
print(f"  nCas9 n={len(nc)} median OR={np.median(nc):.3f} mean={nc.mean():.3f} [{nc.min():.2f},{nc.max():.2f}]")
if len(be)>=3 and len(nc)>=3:
    u,pu=mannwhitneyu(be,nc,alternative='greater')
    print(f"  Mann-Whitney BE4>nCas9 (per-clone OR): U={u:.0f} p={pu:.4f}")
    print(f"  => {'DEAMINASE-SPECIFIC increment (BE4 clones systematically > nCas9)' if pu<0.05 else 'NO clone-level deaminase-specific difference — coupling is SHARED (endogenous APOBEC3), editor adds no separable strand increment'}")
print("  [both cohorts have nCas9/Cas9; only BE4 has deaminase -> this isolates the deaminase, unlike Lei NT-vs-Empty]")
