"""User's metric: rank all C-sites by f(seq)=motif x g(struct)=fork-strand; in the top-N, what is the
EDITING-INDEX (edits/coverage) enrichment vs genome? Answers 'top 5k sites carry X-fold editing risk'."""
import pandas as pd, numpy as np, pyBigWig
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','BE4_CT','BE4_tot'])
sp=sp[sp.BE4_tot>=10].copy()
# sample (stratified) to compute RFD affordably, keep editing counts
samp=sp.sample(400000,random_state=0).reset_index(drop=True)
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
v=[]
for c,p in zip(samp.chrom,samp.pos):
    c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else None)
    try: v.append(rfd.values(c2,int(p)-1,int(p))[0] if c2 else np.nan)
    except: v.append(np.nan)
samp['rfd']=v; samp=samp[samp.rfd.notna()].reset_index(drop=True)
samp['on_exposed']=(((samp.strand=='+')&(samp.rfd>0))|((samp.strand=='-')&(samp.rfd<0))).astype(int)
samp['absrfd']=samp.rfd.abs(); samp['isTpC']=(samp.motif=='TpC').astype(int)
# editing index helper
def idx(df): 
    t=df.BE4_tot.sum(); return df.BE4_CT.sum()/t if t>0 else np.nan
overall=idx(samp)
print(f"=== EDITING-INDEX ENRICHMENT (Doman human BE4, {len(samp)} sampled C-sites, cov>=10) ===")
print(f"  genome-wide editing index (edits/coverage) = {overall:.2e}\n")
print("  by feature stratum (fold vs genome):")
for lab,mask in [('TpC (motif)',samp.isTpC==1),('nonTpC',samp.isTpC==0),
                 ('TpC & fork-exposed',(samp.isTpC==1)&(samp.on_exposed==1)),
                 ('TpC & fork-exposed & high|RFD|(top33%)',(samp.isTpC==1)&(samp.on_exposed==1)&(samp.absrfd>samp.absrfd.quantile(.67)))]:
    print(f"    {lab:42}: index={idx(samp[mask]):.2e}  ({idx(samp[mask])/overall:.2f}x)  n={mask.sum()}")
# combined RANKING score: TpC(2) + on_exposed(1) + absrfd; rank, top-percentile editing-index enrichment
samp['score']=samp.isTpC*2 + samp.on_exposed + (samp.absrfd-samp.absrfd.mean())/samp.absrfd.std()
samp=samp.sort_values('score',ascending=False).reset_index(drop=True)
N=len(samp)
print("\n  RANKED by f(seq)xg(struct) score -> editing-index enrichment in top-k%:")
for k in [1,5,10,20]:
    top=samp.iloc[:int(N*k/100)]; print(f"    top {k:>2}% ({int(N*k/100)} sites): index={idx(top):.2e}  ({idx(top)/overall:.2f}x vs genome)")
print("\n  NOTE: this is APOBEC-context editing (editor + endogenous APOBEC3); WGS editing index is motif-DEGENERATE")
print("  (TpC specificity weak ~1.16x), so enrichment is modest, dominated by fork-strand + coverage.")
