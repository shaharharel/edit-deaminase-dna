"""F1b: full-model vs simple-probe on the RANDOM-negative (0.82-comparable) set.
Feature ablation: seq+-10 / handcrafted-seq / HyenaDNA-256 / all, on the held-out transfer
test (Doman<->McGrath rAPOBEC1). torch-MPS, numpy AUROC."""
import csv, numpy as np, torch, torch.nn as nn
csv.field_size_limit(10**7)
DEV=torch.device("mps" if torch.backends.mps.is_available() else "cpu"); print("device",DEV)
RNG=np.random.default_rng(0)
M={'A':0,'C':1,'G':2,'T':3}; W=256; w10=10
def auroc(y,s):
    y=np.asarray(y); s=np.asarray(s); o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(1,len(s)+1)
    ss=s[o]; i=0
    while i<len(s):
        j=i
        while j+1<len(s) and ss[j+1]==ss[i]: j+=1
        if j>i: r[o[i:j+1]]=(i+1+j+1)/2
        i=j+1
    p=y.sum(); n=len(y)-p
    return float('nan') if p==0 or n==0 else float((r[y==1].sum()-p*(p+1)/2)/(p*n))
# load data
keys=[];src=[];lab=[];seqs=[]
for row in csv.DictReader(open('/tmp/compare_sites.csv')):
    keys.append(row['key']); src.append(row['src']); lab.append(int(row['label'])); seqs.append(row['seq'])
hy=torch.load('/tmp/compare_hyena.pt',map_location='cpu')
keep=[i for i,k in enumerate(keys) if k in hy]
keys=[keys[i] for i in keep]; src=np.array([src[i] for i in keep]); y=np.array([lab[i] for i in keep],float); seqs=[seqs[i] for i in keep]
HYE=np.stack([hy[k].numpy() for k in keys]).astype(np.float32)
print(f"rows {len(keys)} (pos {int(y.sum())}/neg {int((1-y).sum())})")
# seq +-10 one-hot (center=W)
SEQ=np.zeros((len(seqs),(2*w10+1)*4),np.float32)
for i,s in enumerate(seqs):
    sub=s[W-w10:W+w10+1]
    for j,b in enumerate(sub):
        if b in M: SEQ[i,j*4+M[b]]=1
SEQ=SEQ[:,[c for c in range(SEQ.shape[1]) if c//4!=w10]]
# handcrafted-seq: trinuc(16) + GC(3 windows) + kmer3 freq in +-50 (64)
TRI=[a+b for a in 'ACGT' for b in 'ACGT']; TI={t:i for i,t in enumerate(TRI)}
K3=[a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']; KI={k:i for i,k in enumerate(K3)}
HAND=np.zeros((len(seqs),16+3+64),np.float32)
for i,s in enumerate(seqs):
    l,r=s[W-1],s[W+1]
    if l in M and r in M: HAND[i,TI[l+r]]=1
    for wi,win in enumerate([10,50,256]):
        sub=s[W-win:W+win+1]; HAND[i,16+wi]=(sub.count('G')+sub.count('C'))/max(len(sub),1)
    sub=s[W-50:W+51]
    for j in range(len(sub)-2):
        k=sub[j:j+3]
        if k in KI: HAND[i,19+KI[k]]+=1
    HAND[i,19:19+64]/=max(HAND[i,19:19+64].sum(),1)
hyc=torch.load('/tmp/compare_hyena_center.pt',map_location='cpu')
HYC=np.stack([hyc[k].numpy() for k in keys]).astype(np.float32)
FSETS={'seq10':SEQ,'hand-seq':HAND,'hyena-mean':HYE,'hyena-center':HYC,'seq10+hyenaC':np.c_[SEQ,HYC],'all':np.c_[SEQ,HAND,HYE]}
class MLP(nn.Module):
    def __init__(s,d): super().__init__(); s.n=nn.Sequential(nn.Linear(d,128),nn.ReLU(),nn.Dropout(.3),nn.Linear(128,32),nn.ReLU(),nn.Dropout(.3),nn.Linear(32,1))
    def forward(s,x): return s.n(x).squeeze(-1)
def fit_eval(Xtr,ytr,Xte,yte,ep=80):
    mu,sd=Xtr.mean(0),Xtr.std(0)+1e-6; Xtr=(Xtr-mu)/sd; Xte=(Xte-mu)/sd
    Xtr=torch.tensor(Xtr,device=DEV); yt=torch.tensor(ytr,dtype=torch.float32,device=DEV); Xte=torch.tensor(Xte,device=DEV)
    torch.manual_seed(0); m=MLP(Xtr.shape[1]).to(DEV); opt=torch.optim.Adam(m.parameters(),lr=2e-3,weight_decay=1e-4); lf=nn.BCEWithLogitsLoss()
    for _ in range(ep): m.train(); opt.zero_grad(); lf(m(Xtr),yt).backward(); opt.step()
    m.eval()
    with torch.no_grad(): return auroc(yte,torch.sigmoid(m(Xte)).cpu().numpy())
negidx=np.where(y==0)[0]; RNG.shuffle(negidx); h=len(negidx)//2; ntr,nte=negidx[:h],negidx[h:]
D=np.where((src=='Doman_BE4_pilot'))[0]; Mc=np.where((src=='McGrath_2019_iPSC'))[0]
print(f"\n=== F1b feature ablation (random negs, 0.82-comparable). LogReg-probe baseline = 0.82 ===")
print(f"{'featureset':10s} {'dim':>5s} {'Doman->McGr':>12s} {'McGr->Doman':>12s}")
for name,X in FSETS.items():
    md=fit_eval(X[np.r_[D,ntr]],y[np.r_[D,ntr]],X[np.r_[Mc,nte]],y[np.r_[Mc,nte]])
    dm=fit_eval(X[np.r_[Mc,ntr]],y[np.r_[Mc,ntr]],X[np.r_[D,nte]],y[np.r_[D,nte]])
    print(f"{name:10s} {X.shape[1]:5d} {md:12.3f} {dm:12.3f}")
