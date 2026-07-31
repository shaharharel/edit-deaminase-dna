"""McGrath AncBE4max CBE WGS oracle+model panel — the ISOGENIC control (AN21 treated vs ANC uninduced, same line).
Pileup delivered chr1 only (pos->248.9M = chr1 length; chr19/chr22 added 0 rows), so held-out-CHROM is impossible;
use held-out spatial POSITION-BLOCKS within chr1 (train 1st half, rank 2nd half + swap) — valid generalization test.
Schema: chrom,pos,ref,a1t/a1e(AN21_1),a2t/a2e(AN21_2),c1t/c1e(ANC_1),c2t/c2e(ANC_2),tpc. WGS (dilute) CBE, C->T:
ref C -> +strand edited C at center; ref G -> revcomp. Treated=a1+a2, Control=c1+c2. Levanon: NO coverage gate on
the index (Sigma coverage-weighted); light cov>=8 both arms only to define per-site rate for labels/model."""
import pandas as pd, numpy as np, pyfaidx, glob, time
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
t0=time.time(); FLANK=15; BI={'A':0,'C':1,'G':2,'T':3}; comp=np.array([3,2,1,0,4],np.uint8)
d=pd.read_parquet('/mnt/data/mcgrath_wgs/mcgrath_persite.parquet')
d['tt']=d.a1t+d.a2t; d['te']=d.a1e+d.a2e; d['ct']=d.c1t+d.c2t; d['ce']=d.c1e+d.c2e
d=d[(d.tt>=8)&(d.ct>=8)&(d.ref.isin(['C','G']))].reset_index(drop=True)
N=len(d); print(f'McGrath chr1 universe (cov>=8 both) N={N:,} ({time.time()-t0:.0f}s)',flush=True)
te_=d.te.values.astype(float); tt=d.tt.values.astype(float); ce=d.ce.values.astype(float); ct=d.ct.values.astype(float)
trr=te_/np.maximum(tt,1); ctr=ce/np.maximum(ct,1)
ylab=((d.te.values>=2)&(trr>3*ctr+0.002)).astype(int); istpc=d.tpc.values
print(f'edited={int(ylab.sum()):,} ({100*ylab.mean():.3f}%) TpC-frac-of-edits={100*istpc[ylab==1].mean():.0f}%',flush=True)
IDX=lambda ed,tot,i:100*ed[i].sum()/max(tot[i].sum(),1)
def enr(order):
    return '  '.join(f'top{k*100:g}% Tidx={IDX(te_,tt,i):.3f} ratio={IDX(te_,tt,i)/max(IDX(ce,ct,i),1e-9):.1f}x' for k in (0.001,0.01,0.05) for i in [order[:max(int(N*k),1)]])
print(f'(1) ORACLE (rank treated rate): {enr(np.argsort(-trr))}',flush=True)
print(f'    ORACLE editor-specific (rank treated-minus-control): {enr(np.argsort(-(trr-ctr)))}',flush=True)
# ---- subsample to 6M for the model step ONLY (oracle above used full N); keep ALL positives ----
NSS=6_000_000
if N>NSS:
    rng=np.random.default_rng(0)
    pos=np.where(ylab==1)[0]; neg=np.where(ylab==0)[0]
    keepn=rng.choice(neg,size=min(NSS-len(pos),len(neg)),replace=False)
    ss=np.sort(np.concatenate([pos,keepn]))
    d=d.iloc[ss].reset_index(drop=True)
    te_=te_[ss]; tt=tt[ss]; ce=ce[ss]; ct=ct[ss]; trr=trr[ss]; ctr=ctr[ss]; ylab=ylab[ss]; istpc=istpc[ss]
    N=len(d)
    print(f'  subsampled to N={N:,} for model (all {int(ylab.sum())} positives kept; ~67GB->4GB feature matrix)',flush=True)
# seq features (C-oriented), spatial 2-fold by pos median
oh0=0; tri0=120; FDIM=184; X=np.zeros((N,FDIM),np.float32)
fa=pyfaidx.Fasta('/mnt/data/ref/hg38.fa'); lut=np.full(256,4,np.uint8)
for b,i in BI.items(): lut[ord(b)]=i; lut[ord(b.lower())]=i
seq=str(fa['chr1'][:].seq).upper(); codes=lut[np.frombuffer(seq.encode(),np.uint8)]; L=len(codes)
p=d.pos.values.astype(np.int64)-1; v=(p>=FLANK)&(p<L-FLANK)
idx=np.where(v)[0]; p=p[v]; refb=codes[p]; mm=(refb==2)
ctx=np.empty((len(p),2*FLANK+1),np.uint8)
for jj,o in enumerate(range(-FLANK,FLANK+1)): ctx[:,jj]=codes[p+o]
ctx[mm]=comp[ctx[mm][:,::-1]]
flank=np.delete(ctx,FLANK,axis=1)
for j in range(2*FLANK):
    col=flank[:,j]
    for bi in range(4): X[idx[col==bi], oh0+j*4+bi]=1
for i0 in range(ctx.shape[1]-2):
    tri=ctx[:,i0:i0+3]; ok=(tri<4).all(1); c=(tri[:,0]*16+tri[:,1]*4+tri[:,2]).astype(np.int64); np.add.at(X,(idx[ok],tri0+c[ok]),1.0)
X[idx,tri0:tri0+64]/=max(ctx.shape[1]-2,1)
pos_all=d.pos.values; med=np.median(pos_all); fold=(pos_all>med).astype(int)  # spatial held-out halves
print(f'  feats built, spatial fold split at pos {med:.0f} (n0={int((fold==0).sum()):,} n1={int((fold==1).sum()):,}) ({time.time()-t0:.0f}s)',flush=True)
oof=np.full(N,np.nan)
for tef in [0,1]:
    tr=fold!=tef
    m=HistGradientBoostingClassifier(max_iter=400,max_depth=8,learning_rate=0.03,l2_regularization=1.0,class_weight='balanced',random_state=0).fit(X[tr],ylab[tr])
    oof[fold==tef]=m.predict_proba(X[fold==tef])[:,1]
au=roc_auc_score(ylab,oof)
print(f'(2) MODEL GB(flank+tri) held-out-spatial AUROC={au:.4f}: {enr(np.argsort(-oof))}',flush=True)
print('DONE mcgrath_panel2',flush=True)
