"""Control: is the Lei coupling in APOBEC-context (TpC, editor+APOBEC3) or CpG-germline (signature 1)?
Fetch trinuc context from hg19.fa; restrict to TpC; recompute C>T-vs-G>A x sign(RFD) split. Also CpG for contrast."""
import pandas as pd, numpy as np, subprocess
from scipy.stats import fisher_exact
REF="/data/ref/hg19.fa"; BIN="/home/shaharh_quris_ai/miniconda3/envs/apobec/bin"
d=pd.read_parquet('/data/detectseq/lei_editor_sites.parquet')
d=d[d.rfd.notna()&(d.rfd!=0)].copy().reset_index(drop=True)
try:
    import pysam; fa=pysam.FastaFile(REF)
    def ctx(ch,p):  # return + strand bases at p-1,p,p+1 (1-based p)
        try: return fa.fetch(ch,p-2,p+1).upper()
        except: return "NNN"
    print("using pysam",flush=True)
except Exception as e:
    print("pysam unavailable, using samtools faidx batch",flush=True)
    # batch faidx: write regions, call samtools faidx once
    import tempfile,os
    regs="\n".join(f"{c}:{p-1}-{p+1}" for c,p in zip(d.chrom,d.pos))
    r=subprocess.run(f"{BIN}/samtools faidx {REF} -r /dev/stdin",input=regs,shell=True,capture_output=True,text=True)
    seqs={}; cur=None
    for ln in r.stdout.splitlines():
        if ln.startswith(">"): cur=ln[1:].split()[0]
        else: seqs[cur]=seqs.get(cur,"")+ln.upper()
    def ctx(ch,p): return seqs.get(f"{ch}:{p-1}-{p+1}","NNN")
tri=[ctx(c,p) for c,p in zip(d.chrom,d.pos)]
d['tri']=tri
# TpC: for ref C (+strand) the C is middle base -> preceding base (tri[0]) == T
#      for ref G (C on -strand) -> on +strand the base AFTER G (tri[2]) complement is the C's 5' -> tri[2]=='A'
d['isTpC']= ((d.ref=='C')&(d.tri.str[0]=='T')) | ((d.ref=='G')&(d.tri.str[2]=='A'))
d['isCpG']= ((d.ref=='C')&(d.tri.str[2]=='G')) | ((d.ref=='G')&(d.tri.str[0]=='C'))
def splitOR(df):
    ct=df[df.ref=='C']; ga=df[df.ref=='G']
    Cp=int((ct.rfd>0).sum());Cn=int((ct.rfd<0).sum());Gp=int((ga.rfd>0).sum());Gn=int((ga.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,np.nan,len(df)
    o,p=fisher_exact([[Cp,Cn],[Gp,Gn]]); return o,p,len(df)
print(f"total={len(d)}  TpC={d.isTpC.mean():.3f}  CpG={d.isCpG.mean():.3f}  other={1-d.isTpC.mean()-d.isCpG.mean():.3f}",flush=True)
for lab,mask in [('ALL',d.index==d.index),('TpC (APOBEC/editor ctx)',d.isTpC),('CpG (germline sig1 ctx)',d.isCpG),
                 ('non-TpC-non-CpG',~d.isTpC&~d.isCpG)]:
    o,p,n=splitOR(d[mask]); print(f"  {lab:26}: n={n:>8,} OR={o:.3f} p={p:.2g}",flush=True)
d['vaf']=d.mm/d.tot
tp=d[d.isTpC]
print("\n  TpC + VAF-stratified (editor-like low-VAF within APOBEC context):",flush=True)
for lo,hi in [(0,0.10),(0.10,0.20),(0.20,0.35),(0.35,0.65),(0.65,1.01)]:
    o,p,n=splitOR(tp[(tp.vaf>=lo)&(tp.vaf<hi)]); print(f"    VAF[{lo:.2f},{hi:.2f}): n={n:>7,} OR={o:.3f} p={p:.2g}",flush=True)
