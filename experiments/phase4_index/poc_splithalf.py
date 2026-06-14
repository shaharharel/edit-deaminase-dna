"""POC confirmatory test: split-half reproducibility of the per-gene DNA editing index.

The per-clone reproducibility test understates signal (single clones are 4x sparser than the
pool, and correlating single clones on pooled-selected genes is a regression-to-mean trap).
The fair test: pool {clone1,clone5} (half A) vs {clone6,clone7} (half B), compute per-gene
DEI_tpc in each half independently, correlate. If the per-gene index is real and localizable,
the two independent halves should agree.

Reuses the parse logic; re-parses the 4 BE4 clones + Parent once, accumulates per-gene
(edited, tot) counts per clone for tpc and npc, then pools halves.
"""
import os, sys, time
import numpy as np, pandas as pd

DATA='/tmp/poc_dna'; MP=f'{DATA}/mpileups'; BED=f'{DATA}/bed'
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'poc_results')
TREATED=['BE4_clone1','BE4_clone5','BE4_clone6','BE4_clone7']; PARENT='Parent_WGS'
PARENT_VAF_CAP=0.05; MIN_POS=20
def log(*a): print(f'[{time.strftime("%H:%M:%S")}]',*a,file=sys.stderr,flush=True)

log('loading BED...')
pos_index={}; genes=[]; motifs=[]; strands=[]
for bedp,mot in [(f'{BED}/cds_tpc.bed','tpc'),(f'{BED}/cds_npc.bed','npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f=ln.rstrip('\n').split('\t')
            if len(f)<6: continue
            key=(f[0],int(f[1])+1)
            if key in pos_index: continue
            pos_index[key]=len(genes); genes.append(f[3]); motifs.append(mot); strands.append(f[5])
N=len(genes)
genes=np.array(genes,dtype=object); is_tpc=(np.array(motifs)=='tpc'); strand_plus=(np.array(strands)=='+')
log(f'positions {N:,}')

def parse_counts(bases,ref):
    a=c=g=t=0; R=ref.upper(); i=0; n=len(bases)
    while i<n:
        ch=bases[i]
        if ch=='^': i+=2; continue
        if ch=='$': i+=1; continue
        if ch in '+-':
            j=i+1
            while j<n and bases[j].isdigit(): j+=1
            try: i=j+int(bases[i+1:j])
            except ValueError: i+=1
            continue
        if ch=='*': i+=1; continue
        if ch in '.,':
            if R=='A': a+=1
            elif R=='C': c+=1
            elif R=='G': g+=1
            elif R=='T': t+=1
        else:
            u=ch.upper()
            if u=='A': a+=1
            elif u=='C': c+=1
            elif u=='G': g+=1
            elif u=='T': t+=1
        i+=1
    return a,c,g,t

def parse_sample(label):
    tot=np.zeros(N,dtype=np.int32); ed=np.zeros(N,dtype=np.int32)
    with open(f'{MP}/{label}_samtools.mpileup') as fh:
        for ln in fh:
            f=ln.rstrip('\n').split('\t')
            if len(f)<5: continue
            idx=pos_index.get((f[0],int(f[1])))
            if idx is None: continue
            ref=f[2].upper(); a,c,g,t=parse_counts(f[4],ref); tt=a+c+g+t
            if tt<1: continue
            if strand_plus[idx]:
                if ref!='C': continue
                e=t
            else:
                if ref!='G': continue
                e=a
            tot[idx]=tt; ed[idx]=e
    log(f'  parsed {label}')
    return tot,ed

S={s:parse_sample(s) for s in TREATED}
ptot,ped=parse_sample(PARENT)
parent_vaf=np.divide(ped,ptot,out=np.zeros(N),where=ptot>0)
keep=(ptot>=1)&(parent_vaf<=PARENT_VAF_CAP)
c_rate_site=np.divide(ped,ptot,out=np.zeros(N),where=ptot>0)

def half_dei(clones):
    tot=np.zeros(N,dtype=np.int64); ed=np.zeros(N,dtype=np.int64)
    for cl in clones: tot+=S[cl][0]; ed+=S[cl][1]
    k=keep&(tot>=1)
    df=pd.DataFrame({'gene':genes[k],'tpc':is_tpc[k],'ed':ed[k],'tot':tot[k],
                     'c_ed':ped[k],'c_tot':ptot[k]})
    rows=[]
    for gene,sub in df.groupby('gene'):
        m=sub.tpc.values
        for tag,sel in [('tpc',m),('npc',~m)]:
            tt=sub.tot.values[sel].sum(); te=sub.ed.values[sel].sum()
            ct=sub.c_tot.values[sel].sum(); ce=sub.c_ed.values[sel].sum()
            rows.append((gene,tag,sub[sel].shape[0],te/tt if tt else 0,ce/ct if ct else 0))
    r=pd.DataFrame(rows,columns=['gene','tag','npos','t_rate','c_rate'])
    r['dei']=r.t_rate-r.c_rate
    piv=r.pivot(index='gene',columns='tag',values='dei')
    npos=r[r.tag=='tpc'].set_index('gene').npos
    piv['npos_tpc']=npos
    return piv

log('half A = clone1+clone5')
A=half_dei(['BE4_clone1','BE4_clone5'])
log('half B = clone6+clone7')
B=half_dei(['BE4_clone6','BE4_clone7'])
J=A.join(B,lsuffix='_A',rsuffix='_B',how='inner')
J=J[(J.npos_tpc_A>=MIN_POS)&(J.npos_tpc_B>=MIN_POS)]
from scipy.stats import spearmanr,pearsonr
out=[]
out.append("=== SPLIT-HALF reproducibility (pool clone1+5 vs clone6+7) ===")
out.append(f"genes compared: {len(J):,}")
sr=spearmanr(J.tpc_A,J.tpc_B); out.append(f"DEI_tpc  Spearman={sr.correlation:.3f} (p={sr.pvalue:.1e})  Pearson={pearsonr(J.tpc_A,J.tpc_B)[0]:.3f}")
Jc_A=J.tpc_A-J.npc_A; Jc_B=J.tpc_B-J.npc_B
src=spearmanr(Jc_A,Jc_B); out.append(f"DEI_corr Spearman={src.correlation:.3f} (p={src.pvalue:.1e})")
# top-gene reproducibility: does half A's top-500 enrich half B's top?
topA=set(J.sort_values('tpc_A',ascending=False).head(500).index)
topB=set(J.sort_values('tpc_B',ascending=False).head(500).index)
jac=len(topA&topB)/len(topA|topB)
out.append(f"top-500 overlap Jaccard A vs B: {jac:.3f} (random~{500/len(J):.4f}; enrichment {jac/(500/len(J)):.1f}x)")
# does half-B DEI rise across half-A quintiles?
J['qA']=pd.qcut(J.tpc_A.rank(method='first'),5,labels=False)
out.append("mean half-B DEI_tpc by half-A quintile (monotone rise = reproducible signal):")
out.append(J.groupby('qA').tpc_B.mean().round(7).to_string())
txt='\n'.join(out); print(txt)
open(f'{OUT}/splithalf.txt','w').write(txt+'\n')
log('done')
