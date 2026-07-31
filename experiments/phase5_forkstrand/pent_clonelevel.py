import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',columns=['chrom','pos','motif'])
tpc_mask=(sp.motif=='TpC').values
pent=pd.read_parquet('/tmp/poc_dna/tpc_pent.parquet')
# verify alignment (chrom/pos of spectrum-TpC == pent order)
sptpc=sp[tpc_mask].reset_index(drop=True)
assert len(sptpc)==len(pent), (len(sptpc),len(pent))
assert (sptpc.pos.values[:1000]==pent.pos.values[:1000]).all(), "order mismatch"
isY=(pent.pent=='YTCA').values; isR=(pent.pent=='RTCA').values
z=np.load('/tmp/poc_dna/feat/perclone_full.npz',allow_pickle=True)
def clone_pent_rate(clone):
    ct=z[f'{clone}__ct'][:][tpc_mask]; tot=z[f'{clone}__tot'][:][tpc_mask]
    # coverage floor: keep all cov>=1 (index-level, no exclusionary floor per confound a)
    yR=ct[isY].sum()/max(tot[isY].sum(),1); rR=ct[isR].sum()/max(tot[isR].sum(),1)
    return yR, rR, (yR/rR if rR>0 else np.nan)
print("=== CLONE-LEVEL pentamer (YTCA=A3A-like vs RTCA=A3B-like) editing rate, TpC C>T ===")
res={}
for coh,cl in [('BE4',[f'BE4_clone{i}' for i in range(1,9)]),('nCas9',[f'nCas9_clone{i}' for i in range(1,8)])]:
    ratios=[]
    for c in cl:
        yR,rR,ratio=clone_pent_rate(c); ratios.append(ratio)
    res[coh]=np.array(ratios)
    print(f"  {coh}: YTCA/RTCA ratio per clone mean={np.nanmean(res[coh]):.3f} +/- {np.nanstd(res[coh]):.3f}  (n={len(cl)})")
    print(f"       {np.round(res[coh],3)}")
be=res['BE4']; nc=res['nCas9']
t,p=ttest_ind(be,nc,equal_var=False); u,pu=mannwhitneyu(be,nc,alternative='two-sided')
print(f"\n  EDITOR-SPECIFIC pentamer? BE4 YTCA/RTCA={np.nanmean(be):.3f} vs nCas9={np.nanmean(nc):.3f}")
print(f"    Welch p={p:.3f}  MWU p={pu:.3f}")
print("    => if BE4 ratio ~= nCas9 ratio -> editor shares APOBEC3 pentamer pref (pentamer CANNOT separate/enrich editor)")
print("    => if BE4 >> nCas9 -> editor has distinct pentamer fingerprint (could sharpen the f(seq) arm)")
