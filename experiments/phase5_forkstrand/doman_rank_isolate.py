"""Isolate: how much of the Doman ranker 'concentration' is FORK-STRAND biology vs COVERAGE/detectability artifact?"""
import pandas as pd, numpy as np, pyBigWig
from sklearn.ensemble import HistGradientBoostingClassifier
sp=pd.read_parquet('/tmp/poc_dna/feat/per_site_spectrum_24clone.parquet',
    columns=['chrom','pos','motif','strand','BE4_CT','BE4_tot','Parent_CT','Parent_tot'])
sp=sp[(sp.motif=='TpC')&(sp.BE4_tot>=10)&(sp.Parent_tot>=10)&(sp.Parent_CT==0)].copy()
pos=sp[sp.BE4_CT>=3].copy(); neg_pool=sp[sp.BE4_CT==0].copy()
pos['covq']=pd.qcut(pos.BE4_tot,10,labels=False,duplicates='drop')
neg_pool['covq']=pd.qcut(neg_pool.BE4_tot,10,labels=False,duplicates='drop')
negs=[neg_pool[neg_pool.covq==q].sample(min(3*int((pos.covq==q).sum()),len(neg_pool[neg_pool.covq==q])),random_state=0) for q in pos.covq.dropna().unique()]
neg=pd.concat(negs); pos['label']=1; neg['label']=0
d=pd.concat([pos,neg]).reset_index(drop=True)
rfd=pyBigWig.open('/tmp/poc_dna/hg19_hela_rfd.bw'); CHR=set(rfd.chroms())
v=[]
for c,p in zip(d.chrom,d.pos):
    c2=c if c in CHR else ('chr'+c if 'chr'+c in CHR else None)
    try: v.append(rfd.values(c2,int(p)-1,int(p))[0] if c2 else np.nan)
    except: v.append(np.nan)
d['rfd']=v; d=d[d.rfd.notna()&(d.rfd!=0)].reset_index(drop=True)
d['on_exposed']=(((d.strand=='+')&(d.rfd>0))|((d.strand=='-')&(d.rfd<0))).astype(int); d['absrfd']=d.rfd.abs()
y=d.label.values; chrom=d.chrom.values
def locho(cols):
    X=d[cols].astype('float32').values; oof=np.full(len(y),np.nan)
    for c in np.unique(chrom):
        tr=chrom!=c;te=chrom==c
        if y[tr].sum()==0 or y[tr].sum()==tr.sum(): continue
        m=HistGradientBoostingClassifier(max_depth=3,learning_rate=0.05,max_iter=200,min_samples_leaf=30,random_state=0)
        m.fit(X[tr],y[tr]); oof[te]=m.predict_proba(X[te])[:,1]
    return oof
def conc(oof,lab):
    dd=pd.DataFrame({'r':oof,'y':y}).dropna().sort_values('r',ascending=False).reset_index(drop=True)
    N=len(dd);P=dd.y.sum()
    r10=dd.y.iloc[:int(N*.1)].sum()/P; r20=dd.y.iloc[:int(N*.2)].sum()/P
    print(f"  {lab:28}: recall@top10%={r10*100:.1f}% (lift {r10/.1:.2f}x)  @top20%={r20*100:.1f}% (lift {r20/.2:.2f}x)")
print("=== ISOLATING concentration: coverage vs fork-strand (Doman human) ===")
conc(locho(['BE4_tot']),'COVERAGE-only (detectability)')
conc(locho(['on_exposed','absrfd']),'FORK-STRAND-only (callab-immune)')
conc(locho(['BE4_tot','on_exposed','absrfd']),'combined')
print("  => if fork-strand-only lift ~1.0x -> the modest concentration is mostly DETECTABILITY, not editing biology")
