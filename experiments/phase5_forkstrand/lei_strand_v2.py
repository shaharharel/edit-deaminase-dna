"""Human strand test (runs after align) — v2, QA-corrected (ml-code-reviewer).
Decisive human validation of replication-fork-strand coupling of guide-independent CBE editing.

FIXES vs v1:
 (1) Correct mpileup base parsing (parse_bases) — v1 count("T") counted ^mapq/indel tokens, biasing C vs G.
 (2) Empty control requires COVERED (cov>=10) AND unedited (mm==0) at the candidate — not merely 'absent from dict'
     (which was satisfied by Empty being UNCOVERED = a replication-timing/callability confound).
 (3) PRIMARY VERDICT = C>T-vs-G>A 2x2 split (channel x sign(rfd)) Fisher test. This is sign-convention-invariant
     and callability-immune: callability biases act on callable-C-ness (symmetric to whether the event is C>T or
     G>A), so only genuine strand-coupling makes the two channels partition OPPOSITELY along RFD.
     on_exposed vs 0.5 (v1) was NOT adequate (callable-C can be RFD-biased) -> demoted to descriptive + given a
     coverage-matched unedited-negative baseline.
"""
import subprocess, numpy as np, pandas as pd, pyBigWig, random
from scipy.stats import binomtest, fisher_exact, mannwhitneyu
random.seed(0); np.random.seed(0)
BIN="/home/shaharh_quris_ai/miniconda3/envs/apobec/bin"; REF="/data/ref/hg19.fa"
NEG_RESERVOIR=300000   # cap unedited-covered negatives (memory-safe reservoir)

def parse_bases(bases, ref):
    out={'A':0,'C':0,'G':0,'T':0}; R=ref.upper(); i=0; n=len(bases)
    while i<n:
        c=bases[i]
        if c=='^': i+=2; continue
        if c=='$': i+=1; continue
        if c in '+-':
            j=i+1
            while j<n and bases[j].isdigit(): j+=1
            try: i=j+int(bases[i+1:j])
            except ValueError: i+=1
            continue
        if c=='*': i+=1; continue
        if c in '.,':
            if R in out: out[R]+=1
        elif c.upper() in 'ACGT': out[c.upper()]+=1
        i+=1
    return out

def stream_nt(bam):
    """NT: return candidates dict{(c,p):(ref,mm,tot)} for mm>=3, and reservoir of unedited cov>=10 negatives."""
    cand={}; neg=[]; seen_neg=0
    cmd=f"{BIN}/samtools mpileup -f {REF} -q20 -Q20 -d 100000 {bam} 2>/dev/null"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True,bufsize=1)
    for ln in p.stdout:
        f=ln.split("\t")
        if len(f)<5: continue
        ref=f[2].upper()
        if ref not in "CG": continue
        tot=int(f[3])
        if tot<10: continue
        ac=parse_bases(f[4],ref)
        mm=ac['T'] if ref=='C' else ac['A']     # C->T (+strand) or G->A (-strand = C->T on template)
        if mm>=3:
            cand[(f[0],int(f[1]))]=(ref,mm,tot)
        elif mm==0:
            seen_neg+=1                          # reservoir-sample unedited covered sites
            if len(neg)<NEG_RESERVOIR: neg.append((f[0],int(f[1]),ref,tot))
            else:
                k=random.randint(0,seen_neg-1)
                if k<NEG_RESERVOIR: neg[k]=(f[0],int(f[1]),ref,tot)
    return cand,neg

def stream_empty_status(bam, cand):
    """Empty: for candidate positions, record (cov, mm). Memory-safe (only candidate keys)."""
    st={}
    cmd=f"{BIN}/samtools mpileup -f {REF} -q20 -Q20 -d 100000 {bam} 2>/dev/null"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True,bufsize=1)
    for ln in p.stdout:
        f=ln.split("\t")
        if len(f)<5: continue
        key=(f[0],int(f[1]))
        if key not in cand: continue
        ref=f[2].upper(); tot=int(f[3])
        ac=parse_bases(f[4],ref)
        mm=ac['T'] if ref=='C' else (ac['A'] if ref=='G' else 0)
        st[key]=(tot,mm)
    return st

print("streaming NT (candidates + negatives)...",flush=True)
cand,neg=stream_nt("/data/detectseq/SRR11855272.bam")
print(f"  NT candidates(mm>=3)={len(cand):,}  neg reservoir={len(neg):,}",flush=True)
print("streaming Empty (status at candidates)...",flush=True)
emp=stream_empty_status("/data/detectseq/SRR11855254.bam",cand)

# POSITIVES: candidate mm>=3 in NT AND Empty COVERED (cov>=10) & UNEDITED (mm==0)
pos=[(c,p,r,mm,tot) for (c,p),(r,mm,tot) in cand.items()
     if emp.get((c,p),(0,99))[0]>=10 and emp.get((c,p),(0,99))[1]==0]
d=pd.DataFrame(pos,columns=["chrom","pos","ref","mm","tot"])
print(f"  positives (NT mm>=3 & Empty covered-unedited)={len(d):,}",flush=True)

rfd=pyBigWig.open("/data/ref/tracks/hg19_hela_rfd.bw"); CHR=set(rfd.chroms())
def val(ch,p):
    ch2=ch if ch in CHR else ("chr"+ch if "chr"+ch in CHR else None)
    if ch2 is None: return np.nan
    try:
        v=rfd.values(ch2,int(p)-1,int(p))[0]; return v if v==v else np.nan
    except: return np.nan

if len(d)==0:
    print("  NO positives -> abort (check align/coverage)"); raise SystemExit
d["rfd"]=[val(c,p) for c,p in zip(d.chrom,d.pos)]; d=d.dropna(subset=["rfd"])
d=d[d.rfd!=0]
print(f"  positives with nonzero RFD={len(d):,}")

# ---- PRIMARY: C>T-vs-G>A 2x2 split (channel x sign(rfd)) — callability-immune, sign-invariant ----
ct=d[d.ref=="C"]; ga=d[d.ref=="G"]
Cp=int((ct.rfd>0).sum()); Cn=int((ct.rfd<0).sum())
Gp=int((ga.rfd>0).sum()); Gn=int((ga.rfd<0).sum())
table=[[Cp,Cn],[Gp,Gn]]
orr,fp=fisher_exact(table)
fracC=Cp/max(Cp+Cn,1); fracG=Gp/max(Gp+Gn,1)
print("\n=== PRIMARY TEST: C>T vs G>A strand split (callability-immune) ===")
print(f"  C>T (ref C, n={Cp+Cn}): frac at RFD>0 = {fracC:.3f}")
print(f"  G>A (ref G, n={Gp+Gn}): frac at RFD>0 = {fracG:.3f}")
print(f"  2x2 [[C+ {Cp}, C- {Cn}],[G+ {Gp}, G- {Gn}]]  Fisher OR={orr:.2f} p={fp:.3g}")
print(f"  strand-coupling = channels partition OPPOSITELY along RFD; |fracC-fracG|={abs(fracC-fracG):.3f}")

# ---- DESCRIPTIVE: on_exposed (convention rfd>0 => +strand lagging-template) vs matched-neg ----
d["on_exposed"]=(((d.ref=="C")&(d.rfd>0))|((d.ref=="G")&(d.rfd<0)))
oe=d.on_exposed.mean()
# coverage-matched unedited negatives
ndf=pd.DataFrame(neg,columns=["chrom","pos","ref","tot"])
ndf["rfd"]=[val(c,p) for c,p in zip(ndf.chrom,ndf.pos)]; ndf=ndf.dropna(subset=["rfd"]); ndf=ndf[ndf.rfd!=0]
ndf["on_exposed"]=(((ndf.ref=="C")&(ndf.rfd>0))|((ndf.ref=="G")&(ndf.rfd<0)))
# match on coverage decile
d["covq"]=pd.qcut(d.tot,10,labels=False,duplicates="drop")
ndf["covq"]=pd.qcut(ndf.tot,10,labels=False,duplicates="drop")
matched=[]
for q in d.covq.dropna().unique():
    npool=ndf[ndf.covq==q]; k=int((d.covq==q).sum())
    if len(npool)>0: matched.append(npool.sample(min(k,len(npool)),replace=False,random_state=0))
neg_m=pd.concat(matched) if matched else ndf
oe_neg=neg_m.on_exposed.mean()
pval=binomtest(int(d.on_exposed.sum()),len(d),oe_neg if 0<oe_neg<1 else 0.5).pvalue
print("\n=== DESCRIPTIVE: on_exposed vs coverage-matched unedited negatives ===")
print(f"  positives on_exposed={oe:.3f} (n={len(d)}) ; matched-neg baseline={oe_neg:.3f} (n={len(neg_m)})")
print(f"  binom(pos vs neg-baseline) p={pval:.3g}  [note: convention-dependent; PRIMARY split is the real test]")

# ---- VERDICT ----
powered = (Cp+Cn)>=50 and (Gp+Gn)>=50
coupling = fp<0.05 and abs(fracC-fracG)>0.03
if not powered:
    verdict="UNDERPOWERED (too few positives per channel; need deeper NT or looser call)"
elif coupling:
    direction="rfd>0=+strand-lagging (matches mouse)" if fracC>fracG else "INVERTED sign convention vs mouse (coupling still real)"
    verdict=f"HUMAN REPLICATES — strand-coupling significant (Fisher p={fp:.2g}); {direction}"
else:
    verdict="human NULL/weak — no significant C>T-vs-G>A strand partition (callability-immune test)"
print(f"\nVERDICT: {verdict}")
d.to_parquet("/data/detectseq/lei_editor_sites.parquet")
# raw candidates (pre-Empty-filter) + emp status + rfd -> Empty-survival + transcription control locally (QA pre-reg)
cd=pd.DataFrame([(c,p,r,mm,tot) for (c,p),(r,mm,tot) in cand.items()],columns=["chrom","pos","ref","mm","tot"])
cd["emp_cov"]=[emp.get((c,p),(0,-1))[0] for c,p in zip(cd.chrom,cd.pos)]
cd["emp_mm"]=[emp.get((c,p),(0,-1))[1] for c,p in zip(cd.chrom,cd.pos)]
cd["rfd"]=[val(c,p) for c,p in zip(cd.chrom,cd.pos)]
cd.to_parquet("/data/detectseq/lei_candidates.parquet")
print(f"  saved {len(cd)} raw candidates (pre-Empty) for rigorous post-analysis",flush=True)
print("DONE lei_strand")
