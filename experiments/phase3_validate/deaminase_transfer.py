"""Test user hypothesis: motif is DEAMINASE-specific, not treatment-specific => a model trained across
multiple deaminases (with a deaminase feature) should generalize to a HELD-OUT TREATMENT/SOURCE that
reuses a KNOWN deaminase, but NOT to a novel deaminase (zero-shot). Leave-one-SOURCE-out.
Controls (QA): train-only standardization, positives excluded from bg, random + trinuc-matched bg,
enrichment=recall/K, 3 seeds."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
from collections import Counter
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu"); RNG=np.random.default_rng(0)
M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
bn=pd.read_parquet('/tmp/bins_1mb_v3.parquet')
acc={(r.chrom,int(r.bin)):np.array([r.dnase,r.atac,r.rloop,r.h3k27ac],np.float32) for r in bn.itertuples()}
DEAM={'Doman_BE4_pilot':'rAPOBEC1','McGrath_2019_iPSC':'rAPOBEC1','Lei_2021':'rAPOBEC1','Yu_2020':'eng','Richter_2020':'TadA'}
FAMS=['rAPOBEC1','eng','TadA']; BASE={'rAPOBEC1':'C','eng':'C','TadA':'A'}
def famoh(f): v=np.zeros(len(FAMS),np.float32); v[FAMS.index(f)]=1; return v
def ctxget(chrom,pos,base):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if w[W]==COMP[base]: w=rc(w)
    return w if w[W]==base else None
def onehot(w):
    a=np.zeros((2*W+1)*4,np.float32)
    for j,b in enumerate(w):
        if b in M: a[j*4+M[b]]=1
    return a[[k for k in range(len(a)) if k//4!=W]]
def accraw(chrom,pos): return acc.get((chrom,int(pos//1_000_000)), np.full(4,np.nan,np.float32))
def trinuc(w): return w[W-1]+w[W]+w[W+1]
s=pd.read_parquet('data/processed/canonical_sites.parquet'); s=s[s.chrom.isin(CH)]
POS={}  # src -> list of (chrom,pos,oh,acc,fam,trinuc)
allpos=set()
for src,fam in DEAM.items():
    d=s[s.src==src]; base=BASE[fam]; rows=[]
    for r in d.itertuples():
        w=ctxget(r.chrom,int(r.pos),base)
        if w: rows.append((r.chrom,int(r.pos),onehot(w),accraw(r.chrom,int(r.pos)),fam,trinuc(w))); allpos.add((r.chrom,int(r.pos)))
    POS[src]=rows; print(f"{src:20s} {fam:9s} base={base} pos={len(rows)}")
def make_bg(src, matched, mult=4):
    """base-matched random bg for source; trinuc-matched if matched=True; exclude known positives."""
    rows=POS[src]; fam=DEAM[src]; base=BASE[fam]
    tc=Counter(r[5] for r in rows); want={t:c*mult for t,c in tc.items()} if matched else None
    got=Counter(); out=[]; wt=np.array([CL[c] for c in CH],float); wt/=wt.sum(); need=len(rows)*mult; t=0
    while len(out)<need and t<need*80:
        t+=1; c=RNG.choice(CH,p=wt); p=int(RNG.integers(W+2,CL[c]-W-2))
        if (c,p) in allpos: continue
        w=ctxget(c,p,base)
        if not w: continue
        tn=trinuc(w)
        if matched and got[tn]>=want.get(tn,0): continue
        out.append((c,p,onehot(w),accraw(c,p),fam,tn)); got[tn]+=1
    return out
def feats(rows, use_deam):
    O=np.stack([r[2] for r in rows]); A=np.stack([r[3] for r in rows]); D=np.stack([famoh(r[4]) for r in rows])
    return np.c_[O,A,D] if use_deam else np.c_[O,A]
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,64),nn.ReLU(),nn.Dropout(.3),nn.Linear(64,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def auroc(y,sc):
    o=np.argsort(sc);r=np.empty(len(sc));r[o]=np.arange(1,len(sc)+1);p=y.sum();n=len(y)-p
    return float((r[y==1].sum()-p*(p+1)/2)/(p*n))
def enr(y,sc,K=0.10):
    o=np.argsort(-sc); return float(y[o[:int(K*len(y))]].sum()/y.sum())/K
def fit_eval(Xtr,ytr,Xte,yte,seed):
    # train-only standardization; fill nan with train col-mean
    cm=np.nanmean(Xtr,0); ix=np.where(np.isnan(Xtr)); Xtr=Xtr.copy(); Xtr[ix]=np.take(cm,ix[1])
    ixe=np.where(np.isnan(Xte)); Xte=Xte.copy(); Xte[ixe]=np.take(cm,ixe[1])
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6; Xtr=(Xtr-mu)/sd; Xte=(Xte-mu)/sd
    torch.manual_seed(seed); m=MLP(Xtr.shape[1]).to(DEV); opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4)
    pw=torch.tensor([(ytr==0).sum()/max((ytr==1).sum(),1)],dtype=torch.float32,device=DEV); lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    Xs=torch.tensor(Xtr,dtype=torch.float32,device=DEV); ys=torch.tensor(ytr,dtype=torch.float32,device=DEV); idx=np.arange(len(Xs)); rng=np.random.default_rng(seed)
    for e in range(35):
        rng.shuffle(idx)
        for i in range(0,len(idx),4096):
            bi=idx[i:i+4096]; m.train(); opt.zero_grad(); lf(m(Xs[bi]),ys[bi]).backward(); opt.step()
    m.eval()
    with torch.no_grad(): sc=torch.sigmoid(m(torch.tensor(Xte,dtype=torch.float32,device=DEV))).cpu().numpy()
    return auroc(yte,sc),enr(yte,sc)
print("\n=== Leave-one-SOURCE-out: train on OTHER sources (w/ deaminase feature) -> predict held-out treatment ===")
print("    deaminase-in-training? = is the held-out source's deaminase present via a sibling source\n")
BG={src:{'rand':make_bg(src,False),'tri':make_bg(src,True)} for src in DEAM}
order=['Doman_BE4_pilot','McGrath_2019_iPSC','Lei_2021','Yu_2020','Richter_2020']
print(f"{'held-out src':20s} {'deam':9s} {'in-train?':9s} {'AUROC(rand)':>11s} {'enr@10(rand)':>12s} {'AUROC(tri)':>11s} {'enr@10(tri)':>11s}")
for S in order:
    trsrc=[x for x in DEAM if x!=S]
    intrain = DEAM[S] in [DEAM[x] for x in trsrc]
    tr_pos=[r for x in trsrc for r in POS[x]]; tr_bg=[r for x in trsrc for r in BG[x]['rand']]
    Xtr=feats(tr_pos+tr_bg,True); ytr=np.r_[np.ones(len(tr_pos)),np.zeros(len(tr_bg))]
    res={}
    for bgk in ['rand','tri']:
        te=POS[S]+BG[S][bgk]; Xte=feats(te,True); yte=np.r_[np.ones(len(POS[S])),np.zeros(len(BG[S][bgk]))]
        a=[];e=[]
        for seed in (0,1,2):
            au,en=fit_eval(Xtr,ytr,Xte,yte,seed); a.append(au); e.append(en)
        res[bgk]=(np.mean(a),np.mean(e))
    print(f"{S:20s} {DEAM[S]:9s} {('YES' if intrain else 'NO'):9s} {res['rand'][0]:11.3f} {res['rand'][1]:11.1f}x {res['tri'][0]:11.3f} {res['tri'][1]:10.1f}x")
