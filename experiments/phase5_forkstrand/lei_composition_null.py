"""COMPOSITION-NULL control for the Lei 2.46x dose-response (QA scientific-analyst kill-or-firm test).

Confound (agent): high-|RFD| initiation zones are GC/composition-skewed, so the *reference* TpC-C sites and
GpA-G sites (the two channels' denominators) may be asymmetrically distributed w.r.t. RFD sign -- and that
asymmetry can GROW with |RFD|, fabricating the 1.18->1.45->2.46 cross-channel OR gradient with NO strand-coupled
deamination. The C>T-vs-G>A split is callability-immune but NOT composition-immune.

TEST: compute the identical [[C@RFD>0,C@RFD<0],[G@RFD>0,G@RFD<0]] OR on the GENOMIC REFERENCE pool of unedited
TpC-context C sites vs GpA-context G sites (G on + strand followed by A == minus-strand TpC), per |RFD| tertile.
  - If reference OR rises ~1.18->1.45->2.46  => dose-response is COMPOSITION -> 2.46x claim DEAD.
  - If reference OR flat ~1.0 across tertiles => the gradient is real strand-coupled deamination -> claim FIRMS.

No BAM needed: pure hg19.fa composition x RFD. Windowed + subsampled (p) to stay memory-light while rep2 runs.
"""
import pyfaidx, pyBigWig, numpy as np
np.random.seed(0)
FA="/data/ref/hg19.fa"; BW="/data/ref/tracks/hg19_hela_rfd.bw"
CHROMS=["chr1","chr2","chr3","chr7","chr12","chr17"]  # representative; composition skew is genome-general
WIN=5_000_000; P=0.08                                  # subsample prob per qualifying site
fa=pyfaidx.Fasta(FA); bw=pyBigWig.open(BW)
bwchr=set(bw.chroms())
c_rfd=[]; g_rfd=[]                                      # collected RFD at reference TpC-C / GpA-G sites
for ch in CHROMS:
    ch_bw = ch if ch in bwchr else ch.replace("chr","")
    if ch_bw not in bwchr:
        print(f"  {ch}: not in RFD track, skip",flush=True); continue
    L=min(len(fa[ch]), bw.chroms()[ch_bw])
    for s in range(0,L,WIN):
        e=min(s+WIN,L)
        seq=fa[ch][s:e].seq.upper()
        if len(seq)<3: continue
        arr=np.frombuffer(seq.encode(),dtype=np.uint8)          # A65 C67 G71 T84
        rfd=np.array(bw.values(ch_bw,s,e),dtype=np.float32)
        # TpC-C: base==C(67) & prev==T(84).  index i in [1,len-1]
        base=arr[1:-1]; prev=arr[:-2]; nxt=arr[2:]
        r=rfd[1:-1]
        good=np.isfinite(r)&(r!=0)
        tpc = good & (base==67) & (prev==84)
        gpa = good & (base==71) & (nxt==65)             # G followed by A == minus-strand TpC
        # subsample
        ci=np.where(tpc)[0]; gi=np.where(gpa)[0]
        if len(ci): ci=ci[np.random.rand(len(ci))<P]; c_rfd.append(r[ci])
        if len(gi): gi=gi[np.random.rand(len(gi))<P]; g_rfd.append(r[gi])
    print(f"  {ch} done: cumC={sum(len(x) for x in c_rfd):,} cumG={sum(len(x) for x in g_rfd):,}",flush=True)
C=np.concatenate(c_rfd); G=np.concatenate(g_rfd)
def OR(cr,gr):
    Cp=int((cr>0).sum());Cn=int((cr<0).sum());Gp=int((gr>0).sum());Gn=int((gr<0).sum())
    return (Cp*Gn)/max(Cn*Gp,1), (Cp,Cn,Gp,Gn)
allabs=np.abs(np.concatenate([C,G]))
q1,q2=np.quantile(allabs,[1/3,2/3])
print(f"\n=== COMPOSITION-NULL: reference TpC-C vs GpA-G strand OR (NO editing) ===")
print(f"  total ref sites: TpC-C={len(C):,}  GpA-G={len(G):,}  |RFD| tertile cuts=({q1:.2f},{q2:.2f})")
orall,cnt=OR(C,G); print(f"  ALL |RFD|: composition-OR={orall:.3f}  [C+ {cnt[0]}, C- {cnt[1]}, G+ {cnt[2]}, G- {cnt[3]}]")
for lab,lo,hi in [("low",0,q1),("mid",q1,q2),("high",q2,1e9)]:
    cm=(np.abs(C)>=lo)&(np.abs(C)<hi); gm=(np.abs(G)>=lo)&(np.abs(G)<hi)
    o,cnt=OR(C[cm],G[gm]); print(f"  |RFD| {lab:4}: composition-OR={o:.3f} (nC={cm.sum():,} nG={gm.sum():,})")
print("\nINTERPRET: if composition-OR rises ~1.2->1.5->2.5 => the 2.46x dose-response is REFERENCE COMPOSITION, not deamination.")
print("           if composition-OR ~1.0 flat => the edited-site gradient is real strand-coupling.")
