"""Analyze Lei ascertainment-matched composition null: covered-unedited TpC-C vs GpA-G reference sites in the
Detect-seq assay universe. Compute strand OR by |RFD| tertile (FIXED cuts matching edited analysis) -> the correct
null to correct Lei's edited OR. If ~1.0 flat (like Doman genic), Lei edited signal is real; if it rises with |RFD|,
Lei was composition-inflated (genome-wide null was right for Lei's enrichment universe)."""
import pandas as pd, numpy as np, pyfaidx, pyBigWig
FA="/data/ref/hg19.fa"; BW="/data/ref/tracks/hg19_hela_rfd.bw"
fa=pyfaidx.Fasta(FA); bw=pyBigWig.open(BW); bwchr=set(bw.chroms())
d=pd.read_csv("/data/detectseq/covnull/all.tsv",sep="\t",names=["chrom","pos","ref"])
print(f"covered-unedited C/G sites: {len(d):,}")
# TpC context: ref C & prev==T  OR  ref G & next==A (== minus-strand TpC)
def istpc(c,p,ref):
    try:
        if ref=="C": return fa[c][p-2:p-1].seq.upper()=="T"   # prev base (1-based pos)
        else: return fa[c][p:p+1].seq.upper()=="A"            # next base
    except: return False
d["isTpC"]=[istpc(c,p,r) for c,p,r in zip(d.chrom,d.pos,d.ref)]
d=d[d.isTpC].copy()
def rval(c,p):
    cb=c if c in bwchr else c.replace("chr","")
    if cb not in bwchr: return np.nan
    try:
        v=bw.values(cb,int(p)-1,int(p))[0]; return v if v==v else np.nan
    except: return np.nan
d["rfd"]=[rval(c,p) for c,p in zip(d.chrom,d.pos)]
d=d.dropna(subset=["rfd"]); d=d[d.rfd!=0].copy(); d["a"]=d.rfd.abs()
def OR(s):
    cp=s[s.ref=="C"];gp=s[s.ref=="G"]
    Cp=int((cp.rfd>0).sum());Cn=int((cp.rfd<0).sum());Gp=int((gp.rfd>0).sum());Gn=int((gp.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<5: return np.nan,(Cp,Cn,Gp,Gn)
    return (Cp*Gn)/max(Cn*Gp,1),(Cp,Cn,Gp,Gn)
print(f"TpC covered-unedited: {len(d):,}")
print("=== LEI ascertainment-matched composition null (Detect-seq covered-unedited TpC) ===")
for lab,lo,hi in [("ALL",0,9),("low",0,0.25),("mid",0.25,0.54),("high",0.54,9)]:
    s=d[(d.a>=lo)&(d.a<hi)]; o,c=OR(s)
    print(f"  {lab:5}: covered-null OR={o:.3f} (n={len(s)}, {c})")
print("  [compare to genome-wide null 1.06/1.17/1.43 and to Doman genic covered-null ~1.02 flat]")
print("  [Lei EDITED was 1.19/1.45/2.49 -> corrected = edited / THIS null]")
