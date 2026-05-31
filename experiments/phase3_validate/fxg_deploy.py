"""Deployment model: score = f(editor motif spectrum, KNOWN per editor) x g(accessibility, cross-editor).
This is the honest 'generalize between treatments with BE features' model. f from a DISJOINT chrom partition
(non-circular: motif != locations); g from OTHER sources (LOSO). Eval on held-out chrom partition of S.
Shows f-only vs g-only vs fxg, across same-deaminase AND novel-deaminase held-out editors."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from pyfaidx import Fasta
from collections import Counter, defaultdict
HG38="/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"; fa=Fasta(HG38)
W=10; DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu"); RNG=np.random.default_rng(0)
COMP={'A':'T','C':'G','G':'C','T':'A','N':'N'}
def rc(s): return ''.join(COMP[b] for b in reversed(s))
CH=[f'chr{i}' for i in range(1,23)]+['chrX']; CL={c:len(fa[c]) for c in CH}
EVENchr=set(CH[::2])  # f estimated here; eval on ODD
bn=pd.read_parquet('/tmp/bins_1mb_v3.parquet')
acc={(r.chrom,int(r.bin)):np.array([r.dnase,r.atac,r.rloop,r.h3k27ac],np.float32) for r in bn.itertuples()}
amu=np.nanmean(list(acc.values()),0)
DEAM={'Doman_BE4_pilot':'rAPOBEC1','McGrath_2019_iPSC':'rAPOBEC1','Lei_2021':'rAPOBEC1','Yu_2020':'eng','Richter_2020':'TadA'}
BASE={'rAPOBEC1':'C','eng':'C','TadA':'A'}
def ctxget(chrom,pos,base):
    try: w=fa[chrom][pos-1-W:pos+W].seq.upper()
    except: return None
    if len(w)!=2*W+1 or 'N' in w: return None
    if w[W]==COMP[base]: w=rc(w)
    return w if w[W]==base else None
def trinuc(w): return w[W-1]+w[W]+w[W+1]
def accof(chrom,pos):
    a=acc.get((chrom,int(pos//1_000_000)))
    return a if (a is not None and not np.any(np.isnan(a))) else amu
s=pd.read_parquet('data/processed/canonical_sites.parquet'); s=s[s.chrom.isin(CH)]
allpos=set(); SITES=defaultdict(list)
for src,fam in DEAM.items():
    base=BASE[fam]
    for r in s[s.src==src].itertuples():
        w=ctxget(r.chrom,int(r.pos),base)
        if w: SITES[src].append((r.chrom,int(r.pos),trinuc(w),accof(r.chrom,int(r.pos)),base)); allpos.add((r.chrom,int(r.pos)))
def bg_for(base,n):
    out=[];wt=np.array([CL[c] for c in CH],float);wt/=wt.sum();t=0
    while len(out)<n and t<n*60:
        t+=1;c=RNG.choice(CH,p=wt);p=int(RNG.integers(W+2,CL[c]-W-2))
        if (c,p) in allpos: continue
        w=ctxget(c,p,base)
        if w: out.append((c,p,trinuc(w),accof(c,p),base))
    return out
BGc=bg_for('C',60000); BGa=bg_for('A',60000)
def auroc(y,sc):
    o=np.argsort(sc);r=np.empty(len(sc));r[o]=np.arange(1,len(sc)+1);p=y.sum();n=len(y)-p
    return float((r[y==1].sum()-p*(p+1)/2)/(p*n))
def enr(y,sc,K=0.10):
    o=np.argsort(-sc);return float(y[o[:int(K*len(y))]].sum()/y.sum())/K
def spectrum(rows, base):
    """trinuc log-odds vs base composition background (the editor MOTIF; from disjoint chroms)."""
    bg = BGc if base=='C' else BGa
    pc=Counter(r[2] for r in rows); bc=Counter(r[2] for r in bg)
    tot_p=sum(pc.values()); tot_b=sum(bc.values()); lo={}
    for t in set(list(pc)+list(bc)):
        lo[t]=np.log((pc.get(t,0)+1)/(tot_p+1)) - np.log((bc.get(t,0)+1)/(tot_b+1))
    return lo
class GNet(nn.Module):  # g = accessibility-only (deaminase-agnostic)
    def __init__(s):super().__init__();s.n=nn.Sequential(nn.Linear(4,32),nn.ReLU(),nn.Dropout(.3),nn.Linear(32,1))
    def forward(s,x):return s.n(x).squeeze(-1)
def train_g(pos_rows, base):
    bg = BGc if base=='C' else BGa
    X=np.array([r[3] for r in pos_rows]+[r[3] for r in bg],np.float32); y=np.r_[np.ones(len(pos_rows)),np.zeros(len(bg))]
    mu,sd=X.mean(0),X.std(0)+1e-6; Xs=(X-mu)/sd
    torch.manual_seed(0);m=GNet().to(DEV);opt=torch.optim.Adam(m.parameters(),lr=3e-3,weight_decay=1e-4)
    pw=torch.tensor([(y==0).sum()/(y==1).sum()],dtype=torch.float32,device=DEV);lf=nn.BCEWithLogitsLoss(pos_weight=pw)
    Xt=torch.tensor(Xs,device=DEV);yt=torch.tensor(y,dtype=torch.float32,device=DEV)
    for _ in range(250):m.train();opt.zero_grad();lf(m(Xt),yt).backward();opt.step()
    m.eval()
    def score(rows):
        Xr=((np.array([r[3] for r in rows],np.float32))-mu)/sd
        with torch.no_grad(): return torch.sigmoid(m(torch.tensor(Xr,device=DEV))).cpu().numpy()
    return score
print("=== f(known editor motif) x g(cross-editor accessibility): predict HELD-OUT editor's off-targets ===")
print(f"{'held-out editor':20s} {'deam':9s} {'f-only':>7s} {'g-only':>7s} {'fxg':>7s}   (enrichment@top-10%, eval on held-out chroms)")
for S,fam in DEAM.items():
    base=BASE[fam]
    # f: spectrum from EVEN chroms of S (disjoint from eval); g: trained on OTHER sources (LOSO), pooled
    even=[r for r in SITES[S] if r[0] in EVENchr]; odd=[r for r in SITES[S] if r[0] not in EVENchr]
    if len(odd)<100: continue
    spec=spectrum(even, base)
    gpos=[r for x in DEAM if x!=S for r in SITES[x] if BASE[DEAM[x]]==base]  # same-base other sources for g
    if len(gpos)<200: gpos=[r for x in DEAM if x!=S for r in SITES[x]]  # fallback any
    gscore=train_g(gpos, base)
    bg = [r for r in (BGc if base=='C' else BGa) if r[0] not in EVENchr]
    te = odd + bg; y=np.r_[np.ones(len(odd)),np.zeros(len(bg))]
    fsc=np.array([spec.get(r[2], -5.0) for r in te]); gsc=gscore(te)
    fz=(fsc-fsc.mean())/(fsc.std()+1e-6); gz=(gsc-gsc.mean())/(gsc.std()+1e-6)
    print(f"{S:20s} {fam:9s} {enr(y,fz):6.1f}x {enr(y,gz):6.1f}x {enr(y,fz+gz):6.1f}x   (AUROC fxg {auroc(y,fz+gz):.3f}, test pos {len(odd)})")
