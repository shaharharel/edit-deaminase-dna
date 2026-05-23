"""#1 Unified multi-scale net: local conv (f) + regional bin-tracks (g) + editor, trained at SITE level
(more examples), evaluated at BIN level vs hand-split g (pooled 0.47 / Yu 0.64). Clean sources Yu+Lei."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
BIN=1_000_000; W=50; RNG=np.random.default_rng(0); DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
M={'A':0,'C':1,'G':2,'T':3}; COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
bins=pd.read_parquet('/tmp/bins_1mb_v3.parquet').copy(); bins['obs']=bins.ec_Yu+bins.ec_Lei
TR=['dnase','atac','rloop','h3k27ac','mappability','eligible_C','tpc_count']
bins=bins.dropna(subset=TR)
binkey={(r.chrom,int(r.bin)):[getattr(r,t) for t in TR] for r in bins.itertuples()}
tmu=np.array([bins[t].mean() for t in TR]); tsd=np.array([bins[t].std()+1e-6 for t in TR])
def oh(seq):
    X=np.zeros((4,2*W+1),np.float32)
    for i,b in enumerate(seq):
        if b in M: X[M[b],i]=1
    return X
def loc(chrom,pos,flip):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    return rc(w) if flip else w
def make_examples(positions):  # list of (chrom,pos,flip) -> (seq_oh, tracks_norm, chrom, bin)
    S=[];T=[];C=[];B=[]
    for chrom,pos,flip in positions:
        b=pos//BIN; k=(chrom,b)
        if k not in binkey: continue
        w=loc(chrom,pos,flip)
        if w is None: continue
        S.append(oh(w)); T.append((np.array(binkey[k])-tmu)/tsd); C.append(chrom); B.append(b)
    return np.array(S),np.array(T,np.float32),np.array(C),np.array(B)
# positives: Yu+Lei edited C's
sites=pd.read_parquet('data/processed/canonical_sites.parquet')
cl=sites[(sites.edit_class=='CBE')&(sites.src.isin(['Yu_2020','Lei_2021']))]
pos=[(r.chrom,int(r.pos), r.strand=='-') for r in cl.itertuples() if r.chrom in CH]
# random-C negatives (train) + eval random-C set
def rand_C(n):
    out=[];wt=np.array([CL[c] for c in CH],float);wt/=wt.sum();t=0
    while len(out)<n and t<n*30:
        t+=1;c=RNG.choice(CH,p=wt);p=int(RNG.integers(W+2,CL[c]-W-2));bb=fa[c][p-1:p].seq.upper()
        if bb=='C': out.append((c,p,False))
        elif bb=='G': out.append((c,p,True))
    return out
Sp,Tp,Cp,Bp=make_examples(pos)
neg=rand_C(len(pos)); Sn,Tn,Cn,Bn=make_examples(neg)
Sx=np.concatenate([Sp,Sn]);Tx=np.concatenate([Tp,Tn]);Cx=np.concatenate([Cp,Cn]);yx=np.r_[np.ones(len(Sp)),np.zeros(len(Sn))]
ev=rand_C(60000); Se,Te,Ce,Be=make_examples(ev)
print(f"train sites {len(Sx)} (pos {len(Sp)}), eval-C {len(Se)}")
class Net(nn.Module):
    def __init__(s):
        super().__init__()
        s.conv=nn.Sequential(nn.Conv1d(4,32,7,padding=3),nn.ReLU(),nn.Conv1d(32,32,7,padding=3),nn.ReLU(),nn.AdaptiveMaxPool1d(1))
        s.trk=nn.Sequential(nn.Linear(len(TR),16),nn.ReLU())
        s.head=nn.Sequential(nn.Linear(48,32),nn.ReLU(),nn.Dropout(.3),nn.Linear(32,1))
    def forward(s,seq,trk):
        z=s.conv(seq).squeeze(-1); g=s.trk(trk); return s.head(torch.cat([z,g],1)).squeeze(-1)
folds=[CH[i::5] for i in range(5)]  # 5 chromosome folds
binpred={}
for fi,hofold in enumerate(folds):
    tr=~np.isin(Cx,hofold)
    Xs=torch.tensor(Sx[tr],device=DEV);Xt=torch.tensor(Tx[tr],device=DEV);yy=torch.tensor(yx[tr],dtype=torch.float32,device=DEV)
    torch.manual_seed(0);m=Net().to(DEV);opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4);lf=nn.BCEWithLogitsLoss()
    idx=np.arange(len(Xs))
    for ep in range(40):
        RNG.shuffle(idx)
        for i in range(0,len(idx),2048):
            bi=idx[i:i+2048];m.train();opt.zero_grad();lf(m(Xs[bi],Xt[bi]),yy[bi]).backward();opt.step()
    m.eval()
    em=np.isin(Ce,hofold)
    with torch.no_grad():
        pe=torch.sigmoid(m(torch.tensor(Se[em],device=DEV),torch.tensor(Te[em],device=DEV))).cpu().numpy()
    for c,b,pp in zip(Ce[em],Be[em],pe): binpred.setdefault((c,b),[]).append(pp)
bins['pred']=[np.mean(binpred.get((r.chrom,int(r.bin)),[np.nan])) for r in bins.itertuples()]
v=bins.dropna(subset=['pred'])
print(f"\n=== UNIFIED NET (bin-level eval, chrom-holdout) vs hand-split ===")
print(f"  unified net Spearman(pred bin-rate, observed Yu+Lei): {sp(v.pred,v.obs):.3f}")
print(f"  (hand-split g baseline: pooled 0.47, Yu 0.64)")
