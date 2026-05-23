"""Rigorous controls with the WORKING MLP (log1p MSE). Yu+Lei pooled (bulk-only a-priori rule).
delta-over-baseline + block-bootstrap CI + circular-shift null + leave-one-source-out + noise ceiling."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu"); RNG=np.random.default_rng(0)
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values
    return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['obs']=b.ec_Yu+b.ec_Lei; b['logC']=np.log1p(b.eligible_C); b['logTpC']=np.log1p(b.tpc_count)
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def loco(feats,target='obs'):
    pred=np.zeros(len(b))
    for ho in b.chrom.unique():
        tr=b[b.chrom!=ho];te=b[b.chrom==ho]
        mu,sd=tr[feats].mean(0).values,tr[feats].std(0).values+1e-6
        Xtr=((tr[feats].values-mu)/sd).astype(np.float32);Xte=((te[feats].values-mu)/sd).astype(np.float32)
        yt=np.log1p(tr[target].values).astype(np.float32)
        torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
        Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(yt,device=DEV)
        for _ in range(300):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xs),ys).backward();opt.step()
        m.eval()
        with torch.no_grad():pred[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
    return pred
BASE=['logC','logTpC','mappability']; CHROM=['dnase','atac','rloop','h3k27ac']
pb=loco(BASE); pc=loco(BASE+CHROM)
rb,rc=sp(pb,b.obs),sp(pc,b.obs)
print(f"=== Yu+Lei pooled, LOCO Spearman (working MLP) ===")
print(f"  baseline (opportunity+motif+mappability): {rb:.3f}")
print(f"  + chromatin:                              {rc:.3f}   (delta {rc-rb:+.3f})")
chs=b.chrom.unique();boot=[];bd=[]
for _ in range(500):
    pick=RNG.choice(chs,len(chs),replace=True);idx=np.concatenate([np.where(b.chrom.values==c)[0] for c in pick])
    boot.append(sp(pc[idx],b.obs.values[idx]));bd.append(sp(pc[idx],b.obs.values[idx])-sp(pb[idx],b.obs.values[idx]))
print(f"  chromatin 95% CI [{np.percentile(boot,2.5):.3f},{np.percentile(boot,97.5):.3f}]  delta 95% CI [{np.percentile(bd,2.5):+.3f},{np.percentile(bd,97.5):+.3f}]")
null=[sp(pc,np.roll(b.obs.values,int(RNG.integers(50,len(b)-50)))) for _ in range(300)]
print(f"  circular-shift NULL: mean {np.mean(null):+.3f}, 95th {np.percentile(null,95):+.3f}  (real {rc:.3f})")
# leave-one-source-out with MLP: train on ec_Yu, predict, correlate with ec_Lei
pYu=loco(BASE+CHROM,'ec_Yu'); pLei=loco(BASE+CHROM,'ec_Lei')
print(f"  leave-one-source-out: model-of-Yu vs Lei obs: {sp(pYu,b.ec_Lei):.3f} | model-of-Lei vs Yu obs: {sp(pLei,b.ec_Yu):.3f}")
print(f"  [noise ceiling Yu split-half = 0.395; Spearman-Brown full ~0.57]")
print(f"  => chromatin {rc:.3f} vs inter-assay 0.26 and within-Yu-full ~0.57")
