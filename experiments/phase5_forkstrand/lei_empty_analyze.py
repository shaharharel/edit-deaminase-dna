"""EDITOR-ISOLATING control (QA-identified): compute the TpC C>T-vs-G>A x sign(RFD) split OR in the EDITOR-FREE
Empty sample's OWN candidates. Editor claim survives only if OR_TpC(NT=1.61) > OR_TpC(Empty). If Empty~=1.6 -> drop 'editor'."""
import pandas as pd, numpy as np, subprocess, pyBigWig
from scipy.stats import fisher_exact
REF="/data/ref/hg19.fa"; BIN="/home/shaharh_quris_ai/miniconda3/envs/apobec/bin"; D="/data/detectseq"
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
rows=[]
for ln in open(f"{D}/emp_all.txt"):
    f=ln.split("\t")
    if len(f)<5: continue
    ref=f[2].upper()
    if ref not in "CG": continue
    tot=int(f[3])
    if tot<10: continue
    ac=parse_bases(f[4],ref); mm=ac['T'] if ref=='C' else ac['A']
    if mm>=3: rows.append((f[0],int(f[1]),ref,mm,tot))
d=pd.DataFrame(rows,columns=["chrom","pos","ref","mm","tot"])
print(f"EMPTY candidates (mm>=3,cov>=10) = {len(d)}",flush=True)
# RFD
rfd=pyBigWig.open(f"/data/ref/tracks/hg19_hela_rfd.bw"); CHR=set(rfd.chroms())
def val(ch,p):
    ch2=ch if ch in CHR else ("chr"+ch if "chr"+ch in CHR else None)
    if ch2 is None: return np.nan
    try:
        v=rfd.values(ch2,int(p)-1,int(p))[0]; return v if v==v else np.nan
    except: return np.nan
d['rfd']=[val(c,p) for c,p in zip(d.chrom,d.pos)]; d=d[d.rfd.notna()&(d.rfd!=0)].reset_index(drop=True)
# TpC context
regs="\n".join(f"{c}:{p-1}-{p+1}" for c,p in zip(d.chrom,d.pos))
r=subprocess.run(f"{BIN}/samtools faidx {REF} -r /dev/stdin",input=regs,shell=True,capture_output=True,text=True)
seqs={};cur=None
for ln in r.stdout.splitlines():
    if ln.startswith(">"):cur=ln[1:].split()[0]
    else:seqs[cur]=seqs.get(cur,"")+ln.upper()
d['tri']=[seqs.get(f"{c}:{p-1}-{p+1}","NNN") for c,p in zip(d.chrom,d.pos)]
d['isTpC']=((d.ref=='C')&(d.tri.str[0]=='T'))|((d.ref=='G')&(d.tri.str[2]=='A'))
def splitOR(df):
    ct=df[df.ref=='C'];ga=df[df.ref=='G']
    Cp=int((ct.rfd>0).sum());Cn=int((ct.rfd<0).sum());Gp=int((ga.rfd>0).sum());Gn=int((ga.rfd<0).sum())
    if min(Cp+Cn,Gp+Gn)<10: return np.nan,np.nan,len(df),(Cp,Cn,Gp,Gn)
    o,p=fisher_exact([[Cp,Cn],[Gp,Gn]]); return o,p,len(df),(Cp,Cn,Gp,Gn)
oa,pa,na,_=splitOR(d); ot,pt,nt_,cells=splitOR(d[d.isTpC])
print(f"EMPTY (editor-FREE) split: ALL OR={oa:.3f} (n={na}); TpC OR={ot:.3f} p={pt:.2g} (n={nt_}) cells={cells}",flush=True)
print(f"COMPARE: NT_TpC OR=1.609 vs EMPTY_TpC OR={ot:.3f}",flush=True)
print("EDITOR-SPECIFIC iff NT_TpC >> EMPTY_TpC (significant). If EMPTY_TpC~=1.6 -> signal is endogenous APOBEC3, drop 'editor'.",flush=True)
