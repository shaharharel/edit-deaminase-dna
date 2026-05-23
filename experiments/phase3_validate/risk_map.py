"""Steps 2+3: per-editor regional off-target RISK MAP (factored: motif x validated g landscape),
clinical demo for BE4/rAPOBEC1, with cell-type-robustness check + cancer-gene content of top regions.
HONEST: >=1Mb relative-rank screening prior, flag-not-clear, NOT calibrated per-gene."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from bisect import bisect_right
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','dnase_ipsc','dnase_hspc','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['obs']=b.ec_Yu+b.ec_Lei; b['logC']=np.log1p(b.eligible_C); b['logTpC']=np.log1p(b.tpc_count)
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def fit_all(dn):
    feats=['logC','logTpC','mappability',dn,'atac','rloop','h3k27ac']
    mu,sd=b[feats].mean(0).values,b[feats].std(0).values+1e-6
    X=((b[feats].values-mu)/sd).astype(np.float32);y=np.log1p(b.obs.values).astype(np.float32)
    torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
    Xs=torch.tensor(X,device=DEV);ys=torch.tensor(y,device=DEV)
    for _ in range(400):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xs),ys).backward();opt.step()
    m.eval()
    with torch.no_grad():return m(Xs).cpu().numpy()
risk=fit_all('dnase'); b['risk']=risk
# cell-type robustness of the MAP (expert-recommended measured property)
print("=== cell-type robustness of the risk map (rank correlation) ===")
print(f"  HEK293 vs iPSC-DNase map: {sp(risk,fit_all('dnase_ipsc')):.3f}")
print(f"  HEK293 vs HSPC-DNase map: {sp(risk,fit_all('dnase_hspc')):.3f}  (high => map is cell-type-robust)")
# genes per bin + cancer flag (from risk_ranking_v1 which has cancer flag)
rr=pd.read_parquet('data/processed/risk_ranking_v1.parquet')[['gene','cancer']]
canc=set(rr[rr.cancer].gene)
REFGENE="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt"
rg=pd.read_csv(REFGENE,sep='\t',header=None,names=['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','ec','es','ee','sc','name2','a','c','f'])
rg=rg[rg.cdsStart<rg.cdsEnd].copy(); rg['mid']=(rg.txStart+rg.txEnd)//2; rg['bn']=rg.mid//1_000_000
g2bin=rg.groupby(['chrom','bn']).name2.apply(lambda s:sorted(set(s))).to_dict()
b=b.sort_values('risk',ascending=False).reset_index(drop=True)
print("\n=== BE4/rAPOBEC1 TOP-15 regional off-target risk bins (>=1Mb, relative-rank) ===")
print(f"{'region':16s} {'risk':>6s} {'obs':>5s}  cancer-driver genes in region")
for r in b.head(15).itertuples():
    genes=g2bin.get((r.chrom,int(r.bin)),[]); cg=[g for g in genes if g in canc]
    print(f"  {r.chrom+':'+str(int(r.bin))+'Mb':16s} {r.risk:6.2f} {int(r.obs):5d}  {', '.join(cg[:8]) if cg else '(none)'}")
b.to_parquet('data/processed/be4_risk_map.parquet',index=False)
print(f"\nwrote data/processed/be4_risk_map.parquet ({len(b)} bins)")
print("[honest] >=1Mb relative-rank screening prior; flags regions for targeted deep sequencing; cannot CLEAR a region.")
