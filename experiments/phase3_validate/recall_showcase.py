"""SHOWCASE (user reframe): can we PREDICT/RECOVER published guide-independent off-target LOCATIONS?
Train accessibility model -> rank genome by risk -> do HELD-OUT off-target sites fall in top-ranked regions?
Held-out EDITOR (train CBE Yu+Lei -> recover ABE Richter) + held-out SOURCE (train Lei -> recover Yu).
Controls: vs gene-density ranking, vs random. (>=1Mb bins; this is regional recall.)"""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['logC']=np.log1p(b.eligible_C);b['logTpC']=np.log1p(b.tpc_count);b['logA']=np.log1p(1_000_000-b.eligible_C)
rg=pd.read_csv('/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt',sep='\t',header=None,
   names=['bn','nm','chrom','st','ts','te','cs','ce','ec','es','ee','sc','name2','a','c','f'])
rg=rg[rg.cs<rg.ce].copy();rg['bb']=((rg.ts+rg.te)//2)//1_000_000
gd=rg.groupby(['chrom','bb']).name2.nunique().to_dict()
b['gene_density']=[gd.get((x.chrom,int(x.bin)),0) for x in b.itertuples()]
_s=pd.read_parquet('data/processed/canonical_sites.parquet')
def _binct(sub):
    sub=sub.copy();sub['k']=sub.chrom+'_'+(sub.pos//1_000_000).astype(str);c=sub.groupby('k').size().to_dict()
    return np.array([c.get(f"{x.chrom}_{int(x.bin)}",0) for x in b.itertuples()])
b['ec_Richter']=_binct(_s[(_s.edit_class=='ABE')&(_s.src=='Richter_2020')])
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def loco(feats,tgt):
    pred=np.zeros(len(b))
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
    return pred
def recall(rank_vals, site_counts, Ks=(0.05,0.10,0.20,0.30)):
    order=np.argsort(-rank_vals); tot=site_counts.sum(); out={}
    for K in Ks:
        n=int(K*len(rank_vals)); out[K]=site_counts[order[:n]].sum()/tot
    return out
CHROM=['dnase','atac','rloop','h3k27ac']
print("=== SHOWCASE: recover HELD-OUT guide-independent off-target LOCATIONS (regional recall @ top-K% bins) ===")
# 1) train CBE (Yu+Lei) accessibility model -> recover ABE (Richter) sites (held-out editor class)
b['CBE']=b.ec_Yu+b.ec_Lei
risk_cbe=loco(['logC','logTpC','mappability']+CHROM,'CBE')
for name,target,risk in [('ABE Richter (held-out editor-class; model=CBE)', b.ec_Richter.values, risk_cbe),
                         ('CBE Yu (held-out source; model=Lei)', b.ec_Yu.values, loco(['logC','logTpC','mappability']+CHROM,'ec_Lei'))]:
    rmod=recall(risk, target); rgd=recall(b.gene_density.values.astype(float), target)
    print(f"\n  {name}")
    print(f"    {'topK':>6s} {'MODEL recall':>13s} {'gene-dens recall':>17s} {'random':>7s} {'model/random':>13s}")
    for K in (0.05,0.10,0.20,0.30):
        print(f"    {int(K*100):>5d}% {rmod[K]:13.2%} {rgd[K]:17.2%} {K:7.0%} {rmod[K]/K:12.1f}x")
print("\n[model recall > gene-density recall > random => we predict off-target LOCATIONS beyond the confound]")
