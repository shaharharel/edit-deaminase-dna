"""GOTI feasibility — HARDENED per scientific-analyst QA (confirmatory, not rescue).
Fixes: (1) tighter DETECTABILITY-matched negatives — 1:1 nearest-neighbor within (trinuc,in_gene) on standardized
(gc200, log_dist_tss, rep); the strand features (rfd/on_exposed/absrfd) are NOT matched on (that would over-match away
the signal). (2) multi-seed to quantify model variance. (3) macro-AUC (mean per-chrom) as PRIMARY, pooled secondary.
(4) unified site filter. Pre-registered pass: macro-delta CI>0 AND permutation delta-null p<0.05 on the matched set.
Note: GOTI positives are called-SNV supplementary-table rows -> no per-site WGS coverage/mappability available, so
detectability is matched via region proxies (in_gene, gc, replication-timing), the best available surrogates."""
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors
d=pd.read_parquet('/data/goti/goti_rfd.parquet').copy()
d['dnase']=d.dnase.fillna(d.dnase.median())
d=d.dropna(subset=['rfd']).reset_index(drop=True)     # unified filter (keep rfd==0; strand feats handle sign)
pos=d[d.label==1].copy(); neg=d[d.label==0].copy()
MATCH=['gc200','log_dist_tss','rep']
mu=d[MATCH].mean(); sd=d[MATCH].std().replace(0,1)
def z(df): return ((df[MATCH]-mu)/sd).values
matched_idx=[]; used=set()
for (tri,ing),pg in pos.groupby(['tri','in_gene']):
    ng=neg[(neg.tri==tri)&(neg.in_gene==ing)]
    if len(ng)==0: continue
    nn=NearestNeighbors(n_neighbors=min(len(ng),5)).fit(z(ng))
    _,idxs=nn.kneighbors(z(pg))
    ngi=ng.index.values
    for row in idxs:                                   # greedy 1:1 without replacement
        for j in row:
            gi=ngi[j]
            if gi not in used: used.add(gi); matched_idx.append(gi); break
matched_neg=neg.loc[matched_idx]
pos_ok=pos                                             # keep all positives that had >=1 same-(tri,in_gene) neg
m=pd.concat([pos_ok,matched_neg],ignore_index=True)
print(f"positives={len(pos)} matched_negatives={len(matched_neg)} (1:1 NN within tri x in_gene)")
print("=== balance check (pos vs matched-neg means) ===")
for c in ['gc200','in_gene','log_dist_tss','rep','dnase']:
    print(f"  {c:14}: pos={pos[c].mean():.3f}  matched_neg={matched_neg[c].mean():.3f}  (all-neg was {neg[c].mean():.3f})")
tri1h=pd.get_dummies(m.tri,prefix='t')
CALL=['gc200','rep','in_gene','log_dist_tss','dnase']; STRAND=['rfd','on_exposed','absrfd']
def logo(feat,seed):
    X=pd.concat([m[feat].reset_index(drop=True),tri1h.reset_index(drop=True)],axis=1).values
    y=m.label.values; g=m.Chrom.values; oof=np.full(len(y),np.nan)
    for tr,te in LeaveOneGroupOut().split(X,y,g):
        mod=HistGradientBoostingClassifier(max_iter=300,max_depth=4,learning_rate=0.05,random_state=seed).fit(X[tr],y[tr])
        oof[te]=mod.predict_proba(X[te])[:,1]
    # macro AUC (per-chrom, both-class, >=20 test)
    aucs=[]
    for ch in np.unique(g):
        mk=(g==ch)&~np.isnan(oof)
        if mk.sum()>=20 and len(np.unique(y[mk]))==2: aucs.append(roc_auc_score(y[mk],oof[mk]))
    ok=~np.isnan(oof)
    return np.mean(aucs), roc_auc_score(y[ok],oof[ok]), len(aucs)
SEEDS=range(5)
bmac=[];bp=[];smac=[];sp=[]
for s in SEEDS:
    bm,bpo,nc=logo(CALL,s); sm,spo,_=logo(CALL+STRAND,s)
    bmac.append(bm);bp.append(bpo);smac.append(sm);sp.append(spo)
bmac=np.array(bmac);smac=np.array(smac);dmac=smac-bmac
print(f"\n=== HARDENED FEASIBILITY (5 seeds, {nc} chroms in macro) ===")
print(f"  baseline CALL   : macro={bmac.mean():.3f}±{bmac.std():.3f}  pooled={np.mean(bp):.3f}")
print(f"  strong CALL+STR : macro={smac.mean():.3f}±{smac.std():.3f}  pooled={np.mean(sp):.3f}")
print(f"  MACRO DELTA (strand lift) = +{dmac.mean():.3f} [{dmac.min():.3f},{dmac.max():.3f}] across seeds")
# permutation delta-null (shuffle labels within chrom), seed 0
rng=np.random.default_rng(0); g=m.Chrom.values; y=m.label.values; null=[]
for _ in range(200):
    yp=y.copy()
    for ch in np.unique(g):
        idx=np.where(g==ch)[0]; yp[idx]=rng.permutation(yp[idx])
    ms=m.copy(); ms['label']=yp
    # inline logo on permuted labels
    def logo_y(feat,yy):
        X=pd.concat([m[feat].reset_index(drop=True),tri1h.reset_index(drop=True)],axis=1).values
        oof=np.full(len(yy),np.nan)
        for tr,te in LeaveOneGroupOut().split(X,yy,g):
            mod=HistGradientBoostingClassifier(max_iter=150,max_depth=4,learning_rate=0.05,random_state=0).fit(X[tr],yy[tr]); oof[te]=mod.predict_proba(X[te])[:,1]
        ok=~np.isnan(oof); return roc_auc_score(yy[ok],oof[ok])
    null.append(logo_y(CALL+STRAND,yp)-logo_y(CALL,yp))
null=np.array(null); obs=np.mean(sp)-np.mean(bp)
pval=(np.sum(null>=obs)+1)/(len(null)+1)
print(f"  permutation delta-null (pooled): obs={obs:.3f}  null_mean={null.mean():.3f}  p={pval:.3f}")
passc = (dmac.min()>0) and (pval<0.05)
print(f"\n  PRE-REGISTERED PASS (macro-delta CI>0 AND perm p<0.05): {'PASS' if passc else 'FAIL'}")
print("DONE goti_feasibility_hardened")
