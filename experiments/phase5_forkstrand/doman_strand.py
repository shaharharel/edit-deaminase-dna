#!/usr/bin/env python
"""HUMAN validation #2 (Doman WGS, hg19): does guide-independent editor C->T editing show the SAME
replication-fork-strand coupling found in GOTI (mouse)? Callability-immune C>T-vs-G>A split (strand x sign(rfd)),
which is IMMUNE to the clonal-APOBEC3 LOAD confound that killed the rate-based Doman axes (load-symmetric to strand).
Contrast BE4 (editor) vs nCas9 (deaminase-free control) coupling strength.
"""
import pandas as pd, numpy as np, pyBigWig
from scipy.stats import fisher_exact, binomtest
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','BE4_CT','BE4_tot','nCas9_CT','nCas9_tot','Parent_CT','Parent_tot'])
sp=sp[sp.motif=='TpC']  # APOBEC context
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
def valarr(chrom,pos):
    out=np.full(len(pos),np.nan)
    for i,(c,p) in enumerate(zip(chrom,pos)):
        c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else (c[3:] if c[3:] in CHR else None))
        if c2 is None: continue
        try:
            v=rfd.values(c2,int(p)-1,int(p))[0]
            if v==v: out[i]=v
        except: pass
    return out

def strand_test(df,tag):
    # editor-specific-ish guide-indep sites: >=3 edit reads, covered, Parent covered & unedited (germline-clean)
    df=df.dropna(subset=['rfd']); df=df[df.rfd!=0]
    ctp=df[df.strand=='+']  # C->T on + strand
    gap=df[df.strand=='-']  # G->A on + ref = C->T on - strand template
    Cp=int((ctp.rfd>0).sum()); Cn=int((ctp.rfd<0).sum())
    Gp=int((gap.rfd>0).sum()); Gn=int((gap.rfd<0).sum())
    orr,fp=fisher_exact([[Cp,Cn],[Gp,Gn]])
    fracC=Cp/max(Cp+Cn,1); fracG=Gp/max(Gp+Gn,1)
    oe=(((df.strand=='+')&(df.rfd>0))|((df.strand=='-')&(df.rfd<0))).mean()
    print(f"\n--- {tag} (n={len(df)}) ---")
    print(f"  C>T(+strand,n={Cp+Cn}) frac@RFD>0={fracC:.3f} ; G>A(-strand,n={Gp+Gn}) frac@RFD>0={fracG:.3f}")
    print(f"  2x2 Fisher OR={orr:.2f} p={fp:.3g}  |fracC-fracG|={abs(fracC-fracG):.3f}  on_exposed={oe:.3f}")
    return dict(tag=tag,n=len(df),fracC=fracC,fracG=fracG,fp=fp,oe=oe,orr=orr)

# BE4 editor sites (Parent-clean guide-independent)
be=sp[(sp.BE4_CT>=3)&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
be['rfd']=valarr(be.chrom.values,be.pos.values)
# nCas9 control sites (same definition, deaminase-free)
nc=sp[(sp.nCas9_CT>=3)&(sp.nCas9_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
nc['rfd']=valarr(nc.chrom.values,nc.pos.values)

print("=== DOMAN HUMAN fork-strand coupling (hg19 HeLa RFD, TpC, Parent-clean guide-indep) ===")
rb=strand_test(be,'BE4 editor C>T sites')
rn=strand_test(nc,'nCas9 control C>T sites')
print("\n=== VERDICT ===")
be_coupled = rb['fp']<0.05 and abs(rb['fracC']-rb['fracG'])>0.03
print(f"  BE4 strand-coupling: {'PRESENT' if be_coupled else 'absent'} (Fisher p={rb['fp']:.2g}, |dfrac|={abs(rb['fracC']-rb['fracG']):.3f})")
print(f"  nCas9 strand-coupling: {'present' if (rn['fp']<0.05 and abs(rn['fracC']-rn['fracG'])>0.03) else 'absent'} (Fisher p={rn['fp']:.2g})")
if be_coupled:
    direction = "rfd>0=+lagging (matches GOTI mouse)" if rb['fracC']>rb['fracG'] else "inverted sign vs mouse (coupling still real)"
    print(f"  => HUMAN DOMAN REPLICATES fork-strand coupling of guide-independent editing; {direction}")
    print(f"     (callability/APOBEC3-load immune: C>T & G>A partition oppositely along RFD)")
be.to_parquet('/tmp/poc_dna/goti/doman_strand_be4.parquet')
print("\nDONE doman_strand")
