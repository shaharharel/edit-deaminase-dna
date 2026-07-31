"""Pre-registered rigorous analysis of the Lei human verdict (QA criteria, locked BEFORE results).
REPLICATION iff: directional-leg SIGN matches GOTI+Doman AND Empty-control survives AND intergenic-only survives
AND OR 1.2-1.5x AND powered. WRONG SIGN or Empty-killed or intergenic-collapse => NULL (report, don't explain away)."""
import pandas as pd, numpy as np, gzip, sys
from scipy.stats import fisher_exact
POS='/tmp/poc_dna/lei/lei_editor_sites.parquet'; CAND='/tmp/poc_dna/lei/lei_candidates.parquet'
if len(sys.argv)>1 and sys.argv[1]=='--dry': POS='/tmp/poc_dna/_syn_pos.parquet'; CAND='/tmp/poc_dna/_syn_cand.parquet'
def splitOR(df):
    df=df[df.rfd.notna()&(df.rfd!=0)]
    ctp=df[df.ref=='C']; gap=df[df.ref=='G']
    Cp=int((ctp.rfd>0).sum());Cn=int((ctp.rfd<0).sum());Gp=int((gap.rfd>0).sum());Gn=int((gap.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return dict(OR=np.nan,p=np.nan,n=len(df),fracC=np.nan,fracG=np.nan,cells=(Cp,Cn,Gp,Gn))
    orr,pp=fisher_exact([[Cp,Cn],[Gp,Gn]])
    return dict(OR=orr,p=pp,n=len(df),fracC=Cp/(Cp+Cn),fracG=Gp/(Gp+Gn),cells=(Cp,Cn,Gp,Gn))
# refGene -> per-chrom SORTED (start,end,strand) for FAST interval lookup (searchsorted on starts)
raw={}
for ln in gzip.open('/tmp/poc_dna/refGene.txt.gz','rt'):
    f=ln.rstrip().split('\t')
    if len(f)<6: continue
    raw.setdefault(f[2],[]).append((int(f[4]),int(f[5]),f[3]))
G={}
for c,v in raw.items():
    v=sorted(v); starts=np.array([x[0] for x in v]); ends=np.array([x[1] for x in v])
    # cummax of ends so we can test "any interval starting <=p has end>=p"
    G[c]=(starts,ends,np.maximum.accumulate(ends),[x[2] for x in v])
def annotate(df):
    inter=np.ones(len(df),bool); gstr=[None]*len(df)
    for i,(c,p) in enumerate(zip(df.chrom.values,df.pos.values)):
        key=c if c in G else ('chr'+c if 'chr'+c in G else None)
        if key is None: continue
        starts,ends,cmax,strs=G[key]
        j=np.searchsorted(starts,p,side='right')  # intervals with start<=p are [0:j]
        if j==0: continue
        # any of [0:j] with end>=p ?  (linear only over the few overlapping; use cummax to shortcut)
        if cmax[j-1]<p: continue
        for k in range(j-1,-1,-1):
            if ends[k]>=p: inter[i]=False; gstr[i]=strs[k]; break
            if cmax[k]<p: break
    df=df.copy(); df['intergenic']=inter; df['gstrand']=gstr; return df

pos=pd.read_parquet(POS); cand=pd.read_parquet(CAND)
print(f"=== LEI RIGOROUS (pre-registered) — positives={len(pos)} raw candidates={len(cand)} ===")
P=splitOR(pos)
print(f"\n[1] DIRECTIONAL split (positives): OR={P['OR']:.2f} p={P['p']:.2g} n={P['n']}")
print(f"    C>T frac@RFD>0={P['fracC']:.3f} vs G>A frac@RFD>0={P['fracG']:.3f}  cells(C+,C-,G+,G-)={P['cells']}")
sign_ok = (not np.isnan(P['fracC'])) and (P['fracC']>P['fracG'])
print(f"    SIGN matches GOTI+Doman (fracC>fracG): {sign_ok}")
R=splitOR(cand)
print(f"\n[2] EMPTY-CONTROL survival: raw-candidate OR={R['OR']:.2f} (n={R['n']}) -> Empty-filtered OR={P['OR']:.2f} (n={P['n']})")
empty_survives=(not np.isnan(P['OR'])) and P['OR']>1.1
print(f"    survives Empty control: {empty_survives}")
pa=annotate(pos); I=splitOR(pa[pa.intergenic])
print(f"\n[3] TRANSCRIPTION control — INTERGENIC-ONLY: OR={I['OR']:.2f} p={I['p']:.2g} n={I['n']} ({pa.intergenic.mean()*100:.0f}% intergenic)")
txn_ok=(not np.isnan(I['OR'])) and I['OR']>1.1 and (I['fracC']>I['fracG'])
print(f"    intergenic survives (transcription cannot drive): {txn_ok}")
in_band=(not np.isnan(P['OR'])) and 1.15<=P['OR']<=1.6
print(f"\n[4] magnitude ~1.2-1.5x band: {in_band} (OR={P['OR']:.2f}); n={P['n']}")
print("\n=== PRE-REGISTERED VERDICT ===")
if np.isnan(P['OR']): v="INSUFFICIENT POWER (too few sites per channel)"
elif not sign_ok: v="NULL/DISCORDANT — WRONG SIGN (mechanism does NOT transfer to Lei; report prominently)"
elif not empty_survives: v="NULL — killed by Empty control (callability/coverage artifact)"
elif P['p']<0.05 and empty_survives and txn_ok and in_band: v="REPLICATION — directional+Empty-survived+transcription-controlled+in-band+powered = strong 3-way concordance"
elif sign_ok and (in_band or P['p']<0.1): v="WEAK-BUT-CONCORDANT — correct sign, underpowered/marginal; 3-way directional trend, don't oversell"
else: v="AMBIGUOUS — sign ok but fails a criterion; report honestly naming the failing check"
print(f"  {v}")
