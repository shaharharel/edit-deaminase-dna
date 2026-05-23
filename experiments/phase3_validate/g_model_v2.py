"""g v2: (1) does chromatin add beyond opportunity + TpC MOTIF + mappability? (motif-normalized)
(2) cell-type matching: does HEK293 R-loop predict HEK293 sources (Yu/Lei/Chen) better than iPSC (McGrath)?"""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
df=pd.read_parquet('/tmp/bins_1mb_v2.parquet').dropna(subset=['rloop','dnase','h3k27ac','phylop','phastcons','mappability']).reset_index(drop=True)
df['logC']=np.log1p(df.eligible_C); df['logTpC']=np.log1p(df.tpc_count)
print(f"bins {len(df)}, edits {df.edit_count.sum()}")
print("\n=== raw Spearman with pooled edit_count ===")
for f in ['eligible_C','tpc_count','mappability','atac','rloop','dnase','h3k27ac','phylop','phastcons']:
    print(f"  {f:12s} {sp(df[f],df.edit_count):+.3f}")
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def loco(feats,target='edit_count'):
    pred=np.zeros(len(df)); y=np.log1p(df[target].values).astype(np.float32)
    for ho in df.chrom.unique():
        tr=df[df.chrom!=ho]; te=df[df.chrom==ho]
        Xtr=tr[feats].values.astype(np.float32);Xte=te[feats].values.astype(np.float32)
        mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6;Xtr=(Xtr-mu)/sd;Xte=(Xte-mu)/sd
        ytr=np.log1p(tr[target].values).astype(np.float32)
        torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
        Xt=torch.tensor(Xtr,device=DEV);yt=torch.tensor(ytr,device=DEV)
        for _ in range(300):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xt),yt).backward();opt.step()
        m.eval()
        with torch.no_grad():pred[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
    return sp(pred,df[target].values)
print("\n=== TASK 1: does chromatin add beyond opportunity + MOTIF + mappability? (LOCO) ===")
opp=['logC','logTpC','mappability']
chrom=['rloop','dnase','atac','h3k27ac','phylop','phastcons']
print(f"  opportunity+motif+mappability (baseline): {loco(opp):.3f}")
print(f"  + chromatin:                              {loco(opp+chrom):.3f}")
print("\n=== TASK 2: cell-type matching (per-source LOCO, full model w/ HEK293 R-loop) ===")
ctype={'ec_Doman':'HEK293 (orthogonal-assay)','ec_McGrath':'iPSC (MISMATCH)','ec_Yu':'HEK293T','ec_Lei':'HEK293T','ec_Chen':'HEK293T'}
for src,ct in ctype.items():
    if df[src].sum()<50: continue
    r=loco(opp+chrom,target=src)
    print(f"  {src:12s} [{ct:24s}] LOCO Spearman = {r:.3f}  (n_edits={int(df[src].sum())})")
