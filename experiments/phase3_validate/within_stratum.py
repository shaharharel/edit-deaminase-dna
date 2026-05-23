"""Cohort recommendation #1 (control-by-design): does accessibility predict the off-target landscape
WITHIN gene-density strata? If yes, the signal cannot be the gene-density confound."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values
    return float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else np.nan
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['obs']=b.ec_Yu+b.ec_Lei; b['logC']=np.log1p(b.eligible_C); b['logTpC']=np.log1p(b.tpc_count)
rg=pd.read_csv('/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt',sep='\t',header=None,
   names=['bn','name','chrom','strand','ts','te','cs','ce','ec','es','ee','sc','name2','a','c','f'])
rg=rg[rg.cs<rg.ce].copy(); rg['bb']=((rg.ts+rg.te)//2)//1_000_000
gd=rg.groupby(['chrom','bb']).name2.nunique().to_dict()
b['gene_density']=[gd.get((x.chrom,int(x.bin)),0) for x in b.itertuples()]
# chromatin model prediction (LOCO), features = opportunity + chromatin (NO gene_density)
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
feats=['logC','logTpC','mappability','dnase','atac','rloop','h3k27ac']
pred=np.zeros(len(b))
for ho in b.chrom.unique():
    tr=b[b.chrom!=ho];te=b[b.chrom==ho]
    mu,sd=tr[feats].mean(0).values,tr[feats].std(0).values+1e-6
    Xtr=((tr[feats].values-mu)/sd).astype(np.float32);Xte=((te[feats].values-mu)/sd).astype(np.float32)
    yt=np.log1p(tr.obs.values).astype(np.float32)
    torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
    Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(yt,device=DEV)
    for _ in range(300):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xs),ys).backward();opt.step()
    m.eval()
    with torch.no_grad():pred[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
b['pred']=pred
b['gd_q']=pd.qcut(b.gene_density.rank(method='first'),5,labels=False)
print("=== WITHIN gene-density-quintile prediction (control-by-design) ===")
print(f"{'gd quintile':12s} {'n':>5s} {'gene_dens range':>16s} {'pred->obs':>10s} {'DNase->obs':>11s} {'geneDens->obs':>13s}")
for q in range(5):
    s=b[b.gd_q==q]
    print(f"  Q{q+1:<10d} {len(s):5d} {str(int(s.gene_density.min()))+'-'+str(int(s.gene_density.max())):>16s} {sp(s.pred,s.obs):10.3f} {sp(s.dnase,s.obs):11.3f} {sp(s.gene_density,s.obs):13.3f}")
# pooled within-stratum: rank within stratum then correlate
b['pred_w']=b.groupby('gd_q').pred.rank(); b['obs_w']=b.groupby('gd_q').obs.rank()
print(f"\n  POOLED within-stratum Spearman(pred,obs): {sp(b.pred_w,b.obs_w):.3f}")
print(f"  (global pred->obs for reference: {sp(b.pred,b.obs):.3f}; gene_density->obs global: {sp(b.gene_density,b.obs):.3f})")
print("[if within-stratum stays positive across quintiles => accessibility predicts beyond gene density, control-by-design]")
