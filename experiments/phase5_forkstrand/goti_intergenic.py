"""Transcription control for GOTI (mouse), matching the Doman analysis: does the fork-strand signal survive on
INTERGENIC sites (in_gene==0), where transcription-coupled strand asymmetry cannot drive it? Cleanest control."""
import pandas as pd, numpy as np
from scipy.stats import fisher_exact, mannwhitneyu
d=pd.read_parquet('/tmp/poc_dna/goti/goti_rfd.parquet').dropna(subset=['rfd']).copy()
d=d[d.rfd!=0]
print("cols:", [c for c in d.columns]); print("in_gene values:", d.in_gene.value_counts().to_dict())
def report(sub,tag):
    pos=sub[sub.label==1]; neg=sub[sub.label==0]
    oe_p=pos.on_exposed.mean(); oe_n=neg.on_exposed.mean()
    # C>T-vs-G>A x sign(RFD) split on POSITIVES (callability-immune, matches Doman)
    ctp=pos[pos.Cstrand=='+']; gap=pos[pos.Cstrand=='-']
    Cp=int((ctp.rfd>0).sum());Cn=int((ctp.rfd<0).sum());Gp=int((gap.rfd>0).sum());Gn=int((gap.rfd<0).sum())
    orr,fp=fisher_exact([[Cp,Cn],[Gp,Gn]]) if min(Cp+Cn,Gp+Gn)>=10 else (np.nan,np.nan)
    # pos vs neg on_exposed (matched-negative design)
    try: u,pu=mannwhitneyu(pos.on_exposed,neg.on_exposed,alternative='greater')
    except: pu=np.nan
    print(f"\n{tag} (n_pos={len(pos)}, n_neg={len(neg)}):")
    print(f"  on_exposed: pos={oe_p:.3f} vs neg={oe_n:.3f}  (MWU pos>neg p={pu:.2g})")
    print(f"  C>T-vs-G>A split (pos): OR={orr:.2f} p={fp:.2g}  fracC@rfd>0={Cp/max(Cp+Cn,1):.3f} fracG@rfd>0={Gp/max(Gp+Gn,1):.3f}")
report(d,'ALL sites')
report(d[d.in_gene==1],'GENIC (in_gene=1)')
report(d[d.in_gene==0],'INTERGENIC (in_gene=0) — transcription CANNOT drive')
print("\n  => if INTERGENIC retains the signal (pos>neg on_exposed AND split OR>1), GOTI fork-strand is")
print("     replication-attributable (transcription-independent), consistent with the Doman MH transcription control.")
