"""QA-H1: does the chromatin delta survive controlling for GC fraction + GENE DENSITY (not just opportunity)?
If +0.09 collapses with GC+gene-density in the baseline, the 'accessibility' signal is a residual-opportunity confound."""
import numpy as np, pandas as pd, torch, torch.nn as nn
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu")
def sp(a,b):
    a=pd.Series(a).rank().values;b=pd.Series(b).rank().values;return float(np.corrcoef(a,b)[0,1])
b=pd.read_parquet('/tmp/bins_1mb_v3.parquet').dropna(subset=['dnase','atac','rloop','h3k27ac','mappability']).reset_index(drop=True)
b['obs']=b.ec_Yu+b.ec_Lei; b['logC']=np.log1p(b.eligible_C); b['logTpC']=np.log1p(b.tpc_count)
b['gc']=b.eligible_C/1_000_000.0
# gene density per bin (refGene gene midpoints)
rg=pd.read_csv('/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt',sep='\t',header=None,
   names=['bn','name','chrom','strand','ts','te','cs','ce','ec','es','ee','sc','name2','a','c','f'])
rg=rg[rg.cs<rg.ce].copy(); rg['mid']=(rg.ts+rg.te)//2; rg['bb']=rg.mid//1_000_000
gd=rg.groupby(['chrom','bb']).name2.nunique().to_dict()
b['gene_density']=[gd.get((x.chrom,int(x.bin)),0) for x in b.itertuples()]
class MLP(nn.Module):
    def __init__(s,d):super().__init__();s.n=nn.Sequential(nn.Linear(d,32),nn.ReLU(),nn.Dropout(.2),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def loco(feats):
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
    return sp(pred,b.obs.values)
weak=['logC','logTpC','mappability']
strict=weak+['gc','gene_density']
chrom=['dnase','atac','rloop','h3k27ac']
rw=loco(weak); rws=loco(strict)
rcw=loco(weak+chrom); rcs=loco(strict+chrom)
print("=== H1: chromatin delta over WEAK vs STRICT (GC+gene-density) baseline ===")
print(f"  weak baseline (opp+map):              {rw:.3f}")
print(f"  STRICT baseline (+GC+gene-density):   {rws:.3f}")
print(f"  weak + chromatin:                     {rcw:.3f}  (delta over weak  {rcw-rw:+.3f})")
print(f"  STRICT + chromatin:                   {rcs:.3f}  (delta over STRICT {rcs-rws:+.3f})")
print(f"\n  raw Spearman: gene_density vs obs {sp(b.gene_density,b.obs):+.3f}, gc vs obs {sp(b.gc,b.obs):+.3f}, dnase vs obs {sp(b.dnase,b.obs):+.3f}")
print("[if delta-over-STRICT stays ~+0.05-0.09 => accessibility real beyond GC/gene-density; if ->0 => confound]")
