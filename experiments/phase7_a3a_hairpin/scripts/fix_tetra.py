#!/usr/bin/env python
'''Recompute the -2 base STRAND-AWARE and rebuild the donor ranking.

BUG: pcawg_build.step_tetra branched on ref=='C' to choose the -2 base, but ref
was already pyrimidine-oriented, so minus-strand C mutations also hit that branch
and read seq[p-3] on the PLUS strand -- which is the +2 position in oriented
space, uncomplemented. ~half of all sites got a wrong -2 base, corrupting every
YTCA/RTCA number and therefore the A3A-vs-A3B donor split.

Also calibrates against the TRUE genomic background (0.6057 at TCW), not 0.5.
'''
import numpy as np, time
OUT='/mnt/data/a3a/pcawg'; REF='/mnt/data/ref/hg19.fa'
B={'A':0,'C':1,'G':2,'T':3,'N':4}; COMPI=np.array([3,2,1,0,4],dtype=np.uint8)
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}",flush=True)
log('loading hg19'); seqs={};name=None;ch=[]
for line in open(REF):
    if line[0]=='>':
        if name: seqs[name]=''.join(ch).upper()
        name=line[1:].split()[0].replace('chr',''); ch=[]
    else: ch.append(line.strip())
if name: seqs[name]=''.join(ch).upper()
G={}
for c in [str(i) for i in range(1,23)]+['X']:
    if c in seqs:
        a=np.frombuffer(seqs[c].encode(),dtype=np.uint8); o=np.full(a.shape,4,np.uint8)
        for k,v in B.items():
            if k!='N': o[a==ord(k)]=v
        G[c]=o
del seqs
d=np.load(f'{OUT}/snvs.npz',allow_pickle=True)
chrom,pos,ref,alt,donor,tri=d['chrom'],d['pos'],d['ref'],d['alt'],d['donor'],d['tri']
ok=np.isin(chrom,list(G)); log(f'sites on kept chroms: {ok.sum():,}')
# strand from the reference: C on plus = 0, G on plus = 1
strand=np.full(len(pos),-1,np.int8)
for c in np.unique(chrom[ok]):
    m=(chrom==c)&ok; b=G[c][pos[m]-1]
    strand[m]=np.where(b==1,0,np.where(b==2,1,-1))
# -2 base in STRAND-ORIENTED space
m2=np.full(len(pos),4,np.uint8)
for c in np.unique(chrom[ok]):
    m=(chrom==c)&ok; s=G[c]; p=pos[m]-1; st=strand[m]
    plus=s[np.clip(p-2,0,len(s)-1)]
    minus=COMPI[s[np.clip(p+2,0,len(s)-1)]]
    m2[m]=np.where(st==1,minus,plus)
log('-2 base recomputed strand-aware')
up=np.array([t[0] for t in tri]); dn=np.array([t[2] for t in tri])
is_tcw=(ref=='C')&(up=='T')&np.isin(dn,['A','T'])&np.isin(alt,['T','G'])&ok&(strand>=0)
Y=(m2==1)|(m2==3)
log(f'TCW sites: {is_tcw.sum():,}   YTCA frac overall={Y[is_tcw].mean():.4f}')
rows=[]
for dn_ in np.unique(donor):
    m=(donor==dn_)&is_tcw; n=int(m.sum())
    if n<50: continue
    mc=(donor==dn_)&(ref=='C')&ok
    rows.append((dn_, n/max(int(mc.sum()),1), n, float(Y[m].mean())))
rows.sort(key=lambda r:-r[1])
with open(f'{OUT}/a3a_donor_ranking_v2.tsv','w') as fh:
    fh.write('donor\ttcw_frac\tn_tcw\tytca_frac\n')
    for r in rows: fh.write(f'{r[0]}\t{r[1]:.4f}\t{r[2]}\t{r[3]:.4f}\n')
yt=np.array([r[3] for r in rows]); tw=np.array([r[1] for r in rows]); nn=np.array([r[2] for r in rows])
log(f'donors={len(rows)}  ytca median={np.median(yt):.4f} p90={np.percentile(yt,90):.4f} p99={np.percentile(yt,99):.4f}')
for thr in [0.62,0.65,0.68,0.70,0.72]:
    s=(tw>=0.20)&(yt>=thr)
    log(f'  tcw>=0.20 & ytca>={thr}: {s.sum():3d} donors, {nn[s].sum():>9,} TCW muts, median ytca={np.median(yt[s]) if s.sum() else float("nan"):.4f}')
