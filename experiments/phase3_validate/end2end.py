"""End-to-end showcase: consolidated multi-editor accessibility model + CROSS-EDITOR-CLASS generalization.
Train chromatin model on CBE landscape (Yu+Lei) -> predict ABE landscape (Richter), and vice versa.
If it transfers, the accessibility model predicts off-target landscapes editor-class-agnostically."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
s=pd.read_parquet('data/processed/canonical_sites.parquet')
def binct(sub):
    sub=sub.copy(); sub['k']=sub.chrom+'_'+(sub.pos//1_000_000).astype(str); c=sub.groupby('k').size().to_dict()
    return np.array([c.get(f"{x.chrom}_{int(x.bin)}",0) for x in b.itertuples()])
b['CBE']=binct(s[(s.edit_class=='CBE')&(s.src.isin(['Yu_2020','Lei_2021']))])
b['ABE']=binct(s[(s.edit_class=='ABE')&(s.src=='Richter_2020')])
b['logC']=np.log1p(b.eligible_C);b['logA']=np.log1p(1_000_000-b.eligible_C)
CHROM=['dnase','atac','rloop','h3k27ac','mappability']
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def fit_loco(feats,tgt,score_tgt=None):
    """LOCO: train on tgt, predict; correlate with score_tgt (default tgt)."""
    score_tgt=score_tgt or tgt; pred=np.zeros(len(b))
    for ho in b.chrom.unique():
        tr=b[b.chrom!=ho];te=b[b.chrom==ho]
        mu,sd=tr[feats].mean(0).values,tr[feats].std(0).values+1e-6
        Xtr=((tr[feats].values-mu)/sd).astype(np.float32);Xte=((te[feats].values-mu)/sd).astype(np.float32)
        yt=np.log1p(tr[tgt].values).astype(np.float32)
        torch.manual_seed(0);m=MLP(len(feats)).to(DEV);opt=torch.optim.Adam(m.parameters(),lr=5e-3,weight_decay=1e-3)
        Xs=torch.tensor(Xtr,device=DEV);ys=torch.tensor(yt,device=DEV)
        for _ in range(300):m.train();opt.zero_grad();nn.functional.mse_loss(m(Xs),ys).backward();opt.step()
        m.eval()
        with torch.no_grad():pred[te.index]=m(torch.tensor(Xte,device=DEV)).cpu().numpy()
    return sp(pred,b[score_tgt].values)
print("=== consolidated accessibility model: within- and CROSS-editor-class (LOCO Spearman) ===")
print(f"  CBE model -> CBE landscape (within):       {fit_loco(['logC']+CHROM,'CBE'):.3f}")
print(f"  ABE model -> ABE landscape (within):       {fit_loco(['logA']+CHROM,'ABE'):.3f}")
print(f"  CBE model -> ABE landscape (CROSS-CLASS):  {fit_loco(['logC']+CHROM,'CBE','ABE'):.3f}")
print(f"  ABE model -> CBE landscape (CROSS-CLASS):  {fit_loco(['logA']+CHROM,'ABE','CBE'):.3f}")
print(f"\n  raw CBE-landscape vs ABE-landscape correlation: {sp(b.CBE,b.ABE):.3f}")
print("[cross-class transfer ~ within-class => one accessibility model predicts BOTH editors' landscapes]")
