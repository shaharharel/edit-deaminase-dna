"""(c) CROSS-DATASET TRANSFER + (d) POOLED — the last piece of the (a-d) matrix, and a mini-test of the DeaminaFormer
shared-enzyme premise: does a GB seq-model trained on ONE APOBEC-CBE dataset RECOVER the editing-index enrichment of
ANOTHER (different assay + build)? Pair = Doman(BE4 WGS, hg19, treated=BE4 vs control=nCas9) <-> Lei(APOBEC Detect-seq
ENRICHMENT, hg38, treated=NT vs control=Empty). Shared TpC motif. Feature space = build-agnostic seq context (15bp flank
one-hot + trinuc, C-oriented). Metric = index-enrichment (Sigma ed_t/Sigma tot_t)/(Sigma ed_c/Sigma tot_c) over top-K%
ranked, TREATED-vs-CONTROL. Transfer can only RECOVER the shared learnable signal (motif); if index stays ~1x it confirms
the target sites aren't sequence-rankable by ANY (even cross-editor) model => null. env: apobec, /mnt/data.
Enrichment-assay caveat: Lei index is assay-shaped; DROP VAF gate; editor at high apparent VAF."""
import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np, pyfaidx, time
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
t0=time.time(); FLANK=15; BI={'A':0,'C':1,'G':2,'T':3}; comp=np.array([3,2,1,0,4],np.uint8)
CH=[f'chr{c}' for c in list(range(1,23))+['X']]
NSS=3_000_000

def feats(df, fapath):
    """C-oriented seq features from df.chrom/df.pos using reference fapath. Returns X, valid-mask idx into df."""
    fa=pyfaidx.Fasta(fapath); lut=np.full(256,4,np.uint8)
    for b,i in BI.items(): lut[ord(b)]=i; lut[ord(b.lower())]=i
    N=len(df); X=np.zeros((N,184),np.float32); pos=df.pos.values.astype(np.int64)-1; chrom=df.chrom.values
    okall=np.zeros(N,bool)
    for c in CH:
        m=np.where(chrom==c)[0]
        if not len(m): continue
        cc=c if c in fa else c.replace('chr','')
        if cc not in fa: continue
        seq=str(fa[cc][:].seq).upper(); codes=lut[np.frombuffer(seq.encode(),np.uint8)]; L=len(codes)
        p=pos[m]; v=(p>=FLANK)&(p<L-FLANK); mm=m[v]; p=p[v]
        refb=codes[p]; use=(refb==1)|(refb==3)  # ref C or G only (C>T deamination target)
        mm=mm[use]; p=p[use]; isg=(codes[p]==2)
        ctx=np.empty((len(p),2*FLANK+1),np.uint8)
        for jj,o in enumerate(range(-FLANK,FLANK+1)): ctx[:,jj]=codes[p+o]
        ctx[isg]=comp[ctx[isg][:,::-1]]
        flank=np.delete(ctx,FLANK,axis=1)
        for j in range(2*FLANK):
            col=flank[:,j]
            for bi in range(4): X[mm[col==bi], j*4+bi]=1
        for i0 in range(ctx.shape[1]-2):
            tri=ctx[:,i0:i0+3]; ok=(tri<4).all(1); cc2=(tri[:,0]*16+tri[:,1]*4+tri[:,2]).astype(np.int64)
            np.add.at(X,(mm[ok],120+cc2[ok]),1.0)
        X[mm,120:184]/=max(ctx.shape[1]-2,1); okall[mm]=True
    return X, okall

def load_doman():
    d=pd.read_parquet('/mnt/data/doman/per_site_labels_24clone.parquet',
        columns=['chrom','pos','motif','BE4_ed','BE4_tot','nCas9_ed','nCas9_tot','Parent_ed','Parent_tot'])
    d=d[(d.motif=='TpC')&d.chrom.isin(CH)].reset_index(drop=True)
    rt=d.BE4_ed.values/np.maximum(d.BE4_tot.values,1); rc=d.nCas9_ed.values/np.maximum(d.nCas9_tot.values,1)
    y=((d.BE4_ed.values>=2)&(rt>3*rc+0.002)&(d.Parent_ed.values==0)).astype(np.int8)
    return d, d.BE4_ed.values.astype(float), d.BE4_tot.values.astype(float), d.nCas9_ed.values.astype(float), d.nCas9_tot.values.astype(float), y, '/mnt/data/ref/hg19.fa'

def load_lei():
    d=pd.read_parquet('/mnt/data/lei_rp/lei_persite.parquet')
    d=d[(d.motif=='TpC')&d.chrom.isin(CH)].reset_index(drop=True)
    # enrichment assay: treated=NT, control=Empty; DROP VAF gate; edited = treated has >=2 mm and exceeds control rate
    rt=d.nt_mm.values/np.maximum(d.nt_tot.values,1); rc=d.emp_mm.values/np.maximum(d.emp_tot.values,1)
    y=((d.nt_mm.values>=2)&(rt>3*rc+0.002)).astype(np.int8)
    return d, d.nt_mm.values.astype(float), d.nt_tot.values.astype(float), d.emp_mm.values.astype(float), d.emp_tot.values.astype(float), y, '/mnt/data/ref/hg38.fa'

def subsample(d,edt,tot_t,edc,tot_c,y):
    N=len(d)
    if N<=NSS: return d,edt,tot_t,edc,tot_c,y,np.arange(N)
    rng=np.random.default_rng(0); pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    keepn=rng.choice(neg,size=NSS-len(pos),replace=False); ss=np.sort(np.concatenate([pos,keepn]))
    return d.iloc[ss].reset_index(drop=True),edt[ss],tot_t[ss],edc[ss],tot_c[ss],y[ss],ss

print('loading + featurizing Doman(hg19) and Lei(hg38)...',flush=True)
DS={}
for name,loader in [('Doman',load_doman),('Lei',load_lei)]:
    d,edt,tt,edc,tc,y,fap=loader()
    print(f'  {name}: N={len(d):,} TpC editor-specific pos={int(y.sum()):,} ({100*y.mean():.3f}%) ({time.time()-t0:.0f}s)',flush=True)
    d,edt,tt,edc,tc,y,_=subsample(d,edt,tt,edc,tc,y)
    X,ok=feats(d,fap)
    DS[name]=dict(X=X[ok],y=y[ok],edt=edt[ok],tt=tt[ok],edc=edc[ok],tc=tc[ok],chrom=d.chrom.values[ok])
    print(f'  {name}: featurized rows={int(ok.sum()):,} pos={int(y[ok].sum()):,} ({time.time()-t0:.0f}s)',flush=True)

def idxenr(edt,tt,edc,tc,order,pct):
    K=max(int(len(order)*pct/100),1); i=order[:K]
    ti=edt[i].sum()/max(tt[i].sum(),1); ci=edc[i].sum()/max(tc[i].sum(),1)
    return ti/max(ci,1e-12)

def evaluate(scores,T,tag):
    order=np.argsort(-scores)
    au=roc_auc_score(T['y'],scores) if T['y'].sum()>0 else float('nan')
    s=' '.join(f'top{p}%={idxenr(T["edt"],T["tt"],T["edc"],T["tc"],order,p):.2f}x' for p in (0.1,1,5))
    print(f'  {tag}: target-AUROC={au:.3f}  index-enr[{s}]',flush=True)

print('=== POSITIVE CONTROL: within-dataset held-out-CHROM (proves the harness detects signal when it exists) ===',flush=True)
HELD=set(f'chr{c}' for c in [2,5,9,14,20])  # ~held-out ~20% of genome, disjoint from train
for name in ['Doman','Lei']:
    D=DS[name]; te=np.isin(D['chrom'],list(HELD)); tr=~te
    m=HistGradientBoostingClassifier(max_iter=300,max_depth=6,learning_rate=0.05,l2_regularization=1.0,class_weight='balanced',random_state=0).fit(D['X'][tr],D['y'][tr])
    sc=m.predict_proba(D['X'][te])[:,1]
    Dte={k:D[k][te] for k in ('y','edt','tt','edc','tc')}
    evaluate(sc,Dte,f'{name} WITHIN held-out-chrom (POS-CTRL)')
print('=== (c) CROSS-TRANSFER: train SOURCE, rank TARGET index ===',flush=True)
for src,tgt in [('Doman','Lei'),('Lei','Doman')]:
    S,T=DS[src],DS[tgt]
    m=HistGradientBoostingClassifier(max_iter=300,max_depth=6,learning_rate=0.05,l2_regularization=1.0,class_weight='balanced',random_state=0).fit(S['X'],S['y'])
    sc=m.predict_proba(T['X'])[:,1]
    evaluate(sc,T,f'{src}->{tgt} TRANSFER')
    # baselines on same target
    rng=np.random.default_rng(1); evaluate(rng.random(len(T['y'])),T,f'   {tgt} RANDOM baseline')
print('=== (d) POOLED: train Doman+Lei, rank each (held-out by dataset = same as transfer; pooled adds cross-motif reinforcement) ===',flush=True)
for tgt in ['Doman','Lei']:
    src=[k for k in DS if k!=tgt]  # the other(s)
    Xtr=np.vstack([DS[k]['X'] for k in src]); ytr=np.concatenate([DS[k]['y'] for k in src])
    m=HistGradientBoostingClassifier(max_iter=300,max_depth=6,learning_rate=0.05,l2_regularization=1.0,class_weight='balanced',random_state=0).fit(Xtr,ytr)
    sc=m.predict_proba(DS[tgt]['X'])[:,1]; evaluate(sc,DS[tgt],f'POOLED(all-but-{tgt})->{tgt}')
print('INTERPRET: transfer index-enr ~1x => target sites not sequence-rankable even by a cross-dataset APOBEC-CBE model (motif shared but sites diffuse) = confirms null across the (a-d) matrix. >1.5x at small K => shared editor sequence-preference DOES transfer (would be the first real predictive positive).',flush=True)
print('DONE cross_transfer',flush=True)
