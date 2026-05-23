"""Symmetric cell-type test: does the MATCHED cell-type DNase predict each source best?
McGrath=iPSC should prefer iPSC-DNase; Yu/Lei=HEK293 should prefer HEK293-DNase."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
df=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','dnase_ipsc','dnase_hspc','mappability']).reset_index(drop=True)
df['logC']=np.log1p(df.eligible_C);df['logTpC']=np.log1p(df.tpc_count)
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,16),nn.ReLU(),nn.Dropout(.2),nn.Linear(16,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def loco(feats,target):
    pred=np.zeros(len(df))
    for ho in df.chrom.unique():
        tr=df[df.chrom!=ho];te=df[df.chrom==ho]
        Xtr=tr[feats].values.astype(np.float32);Xte=te[feats].values.astype(np.float32)
        mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6;Xtr=(Xtr-mu)/sd;Xte=(Xte-mu)/sd
        yt=np.log1p(tr[target].values).astype(np.float32)
        torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
        Xt=torch.tensor(Xtr,device=DEV);yy=torch.tensor(yt,device=DEV)
        for _ in range(250):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xt),yy).backward();opt.step()
        m.eval()
        with torch.no_grad():pred[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
    return sp(pred,df[target].values)
base=['logC','logTpC','mappability']
print(f"{'source(celltype)':24s} {'+HEK293':>9s} {'+iPSC':>9s} {'+HSPC':>9s}  matched-wins?")
rows=[('ec_McGrath','iPSC'),('ec_Yu','HEK293'),('ec_Lei','HEK293'),('ec_Doman','HEK293-ortho')]
for src,ct in rows:
    if df[src].sum()<50: continue
    hek=loco(base+['dnase'],src); ips=loco(base+['dnase_ipsc'],src); hsp=loco(base+['dnase_hspc'],src)
    best='HEK293' if hek>=max(ips,hsp) else ('iPSC' if ips>=hsp else 'HSPC')
    flag='<-- MATCHED' if (ct.startswith('iPSC') and best=='iPSC') or (ct.startswith('HEK293') and best=='HEK293') else ''
    print(f"  {src+' ('+ct+')':24s} {hek:9.3f} {ips:9.3f} {hsp:9.3f}  best={best} {flag}")
print("\n[if McGrath prefers iPSC-DNase and Yu/Lei prefer HEK293-DNase => cell-type matching confirmed both ways]")
