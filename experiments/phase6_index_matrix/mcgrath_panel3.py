"""McGrath CpG-vs-TpC DECOMPOSITION — is the AUROC=0.664 / top1% 1.4x model predicting EDITOR (TpC) off-targets, or
just GERMLINE C>T (CpG deamination, trivially sequence-predictable)? The edited-call set is only 28% TpC => germline
suspected. DECISIVE checks: (A) motif composition of edited calls (%TpC vs %CpG vs background); (B) motif composition
of the model's TOP-1% ranked panel (CpG-enriched => predicting germline); (C) retrain model on TpC-EDITOR-ONLY positives
-> does AUROC/enrichment survive, or collapse to ~0.5 (=> the 0.664 was entirely CpG/germline). Held-out spatial 2-fold
within chr1 (chr1-only pileup; NOT held-out-chrom => regional-leakage caveat stands). env: apobec, /mnt/data."""
import pandas as pd, numpy as np, pyfaidx, time
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
t0=time.time(); FLANK=15; BI={'A':0,'C':1,'G':2,'T':3}; comp=np.array([3,2,1,0,4],np.uint8)
d=pd.read_parquet('/mnt/data/mcgrath_wgs/mcgrath_persite.parquet')
d['tt']=d.a1t+d.a2t; d['te']=d.a1e+d.a2e; d['ct']=d.c1t+d.c2t; d['ce']=d.c1e+d.c2e
d=d[(d.tt>=8)&(d.ct>=8)&(d.ref.isin(['C','G']))].reset_index(drop=True)
N=len(d); print(f'McGrath chr1 universe (cov>=8 both) N={N:,} ({time.time()-t0:.0f}s)',flush=True)
te_=d.te.values.astype(float); tt=d.tt.values.astype(float); ce=d.ce.values.astype(float); ct=d.ct.values.astype(float)
trr=te_/np.maximum(tt,1); ctr=ce/np.maximum(ct,1)
ylab=((d.te.values>=2)&(trr>3*ctr+0.002)).astype(int)
# ---- subsample 6M keeping all positives ----
NSS=6_000_000
if N>NSS:
    rng=np.random.default_rng(0); pos=np.where(ylab==1)[0]; neg=np.where(ylab==0)[0]
    keepn=rng.choice(neg,size=NSS-len(pos),replace=False); ss=np.sort(np.concatenate([pos,keepn]))
    d=d.iloc[ss].reset_index(drop=True); te_=te_[ss];tt=tt[ss];ce=ce[ss];ct=ct[ss];trr=trr[ss];ctr=ctr[ss];ylab=ylab[ss]; N=len(d)
    print(f'  subsampled N={N:,} (all {int(ylab.sum())} positives kept)',flush=True)
# ---- C-oriented context: m1 (5prime, -1), p1 (3prime, +1) ----
fa=pyfaidx.Fasta('/mnt/data/ref/hg38.fa'); lut=np.full(256,4,np.uint8)
for b,i in BI.items(): lut[ord(b)]=i; lut[ord(b.lower())]=i
seq=str(fa['chr1'][:].seq).upper(); codes=lut[np.frombuffer(seq.encode(),np.uint8)]; L=len(codes)
p=d.pos.values.astype(np.int64)-1; v=(p>=FLANK)&(p<L-FLANK); idx=np.where(v)[0]; p=p[v]; refb=codes[p]; mm=(refb==2)
ctx=np.empty((len(p),2*FLANK+1),np.uint8)
for jj,o in enumerate(range(-FLANK,FLANK+1)): ctx[:,jj]=codes[p+o]
ctx[mm]=comp[ctx[mm][:,::-1]]  # G-ref -> revcomp to C-orientation
m1=np.full(N,4,np.uint8); p1=np.full(N,4,np.uint8); m1[idx]=ctx[:,FLANK-1]; p1[idx]=ctx[:,FLANK+1]
istpc=(m1==3); iscpg=(p1==2)   # TpC editor motif (T at -1) ; CpG germline motif (G at +1)
def comp_str(m):
    n=max(m.sum(),1); return f'TpC={100*istpc[m].mean():.0f}% CpG={100*iscpg[m].mean():.0f}% (n={int(m.sum()):,})'
print('=== (A) MOTIF COMPOSITION ===',flush=True)
print(f'  edited calls     : {comp_str(ylab==1)}',flush=True)
print(f'  background(all C) : {comp_str(ylab==0)}',flush=True)
# ---- features (flank one-hot + trinuc), spatial 2-fold ----
oh0=0; tri0=120; FDIM=184; X=np.zeros((N,FDIM),np.float32)
flank=np.delete(ctx,FLANK,axis=1)
for j in range(2*FLANK):
    col=flank[:,j]
    for bi in range(4): X[idx[col==bi], oh0+j*4+bi]=1
for i0 in range(ctx.shape[1]-2):
    tri=ctx[:,i0:i0+3]; ok=(tri<4).all(1); c=(tri[:,0]*16+tri[:,1]*4+tri[:,2]).astype(np.int64); np.add.at(X,(idx[ok],tri0+c[ok]),1.0)
X[idx,tri0:tri0+64]/=max(ctx.shape[1]-2,1)
med=np.median(d.pos.values); fold=(d.pos.values>med).astype(int)
IDX=lambda ed,tot,i:100*ed[i].sum()/max(tot[i].sum(),1)
def enr(order):
    return '  '.join(f'top{k*100:g}% ratio={IDX(te_,tt,i)/max(IDX(ce,ct,i),1e-9):.1f}x' for k in (0.001,0.01,0.05) for i in [order[:max(int(N*k),1)]])
def run_model(ylab_use,tag):
    oof=np.full(N,np.nan)
    for tef in [0,1]:
        tr=fold!=tef
        if ylab_use[tr].sum()<20: print(f'  {tag}: too few positives in train fold, skip',flush=True); return None
        m=HistGradientBoostingClassifier(max_iter=400,max_depth=8,learning_rate=0.03,l2_regularization=1.0,class_weight='balanced',random_state=0).fit(X[tr],ylab_use[tr])
        oof[fold==tef]=m.predict_proba(X[fold==tef])[:,1]
    au=roc_auc_score(ylab_use,oof); order=np.argsort(-oof)
    top1=order[:max(int(N*0.01),1)]
    print(f'  {tag}: AUROC={au:.4f}  index-enr[{enr(order)}]  TOP1%-panel-motif[TpC={100*istpc[top1].mean():.0f}% CpG={100*iscpg[top1].mean():.0f}%]',flush=True)
    return au
print('=== (B/C) MODEL under different POSITIVE-label definitions (same features, same spatial folds) ===',flush=True)
print('  [if TpC-only AUROC collapses ~0.5 and top-panel is CpG-heavy => model predicts GERMLINE not EDITOR]',flush=True)
run_model(ylab,'all-edited (orig 0.664 repro)')
run_model((ylab==1)&istpc,'TpC-EDITOR-only positives')
run_model((ylab==1)&iscpg,'CpG-GERMLINE-only positives')
print('DONE mcgrath_panel3',flush=True)
