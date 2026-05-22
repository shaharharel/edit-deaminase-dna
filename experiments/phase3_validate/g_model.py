"""g-landscape model: does chromatin predict the >=1Mb guide-indep CBE off-target landscape
BEYOND opportunity (eligible-C) + ascertainment (mappability)? Leave-one-chromosome-out."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def spearman(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values
    return float(np.corrcoef(a,b)[0,1])
df=pd.read_parquet('/tmp/bins_1mb.parquet').dropna(subset=['rloop','dnase','phylop','phastcons','mappability']).reset_index(drop=True)
df['y']=np.log1p(df.edit_count); df['logC']=np.log1p(df.eligible_C)
print(f"bins {len(df)}, total edits {df.edit_count.sum()}")
# raw associations
print("\n=== raw Spearman with edit_count ===")
for f in ['eligible_C','mappability','rloop','dnase','phylop','phastcons']:
    print(f"  {f:12s} {spearman(df[f],df.edit_count):+.3f}")
# label reproducibility across sources at 1Mb (target validity)
srcs=[c for c in df.columns if c.startswith('ec_')]
print("\n=== label cross-source Spearman (target reproducibility) ===")
for i in range(len(srcs)):
    for j in range(i+1,len(srcs)):
        print(f"  {srcs[i]} vs {srcs[j]}: {spearman(df[srcs[i]],df[srcs[j]]):+.3f}")
class MLP(nn.Module):
    def __init__(s,d): super().__init__(); s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x): return s.n(x).squeeze(-1)
def loco(feats):
    preds=np.zeros(len(df)); chroms=df.chrom.unique()
    for ho in chroms:
        tr=df[df.chrom!=ho]; te=df[df.chrom==ho]
        Xtr=tr[feats].values.astype(np.float32); Xte=te[feats].values.astype(np.float32)
        mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6; Xtr=(Xtr-mu)/sd; Xte=(Xte-mu)/sd
        ytr=tr.y.values.astype(np.float32)
        torch.manual_seed(0); m=MLP(len(feats)).to(DEV); opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
        Xt=torch.tensor(Xtr,device=DEV); yt=torch.tensor(ytr,device=DEV)
        for _ in range(300): m.train(); opt.zero_grad(); nn.functional.mse_loss(m(Xt),yt).backward(); opt.step()
        m.eval()
        with torch.no_grad(): preds[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
    return spearman(preds, df.edit_count.values)
print("\n=== leave-one-chromosome-out Spearman(pred, observed edit_count) ===")
base=['logC','mappability']
print(f"  baseline (eligible_C + mappability):        {loco(base):.3f}")
print(f"  + R-loop:                                   {loco(base+['rloop']):.3f}")
print(f"  + chromatin (rloop,dnase,phylop,phastcons): {loco(base+['rloop','dnase','phylop','phastcons']):.3f}")
print("\n[interpretation] if +chromatin >> baseline, the g landscape is real & chromatin-driven.")
print("if ~equal, edit landscape is just opportunity+ascertainment (no extra chromatin signal).")
