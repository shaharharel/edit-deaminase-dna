#!/usr/bin/env python
"""IS THE GENOMIC-SCALE CANCER RESULT IN-SAMPLE? The obvious attack on section 34.

s14 fits ONE model on ALL of a3a_trainset_v5.npz and then scores the whole genome, joining
PCAWG mutations FROM THE SAME 21 DONORS the model trained on. The trainset holds 83,999
positives drawn from those donors; the genomic join found 92,823 mutations from them.

IF THOSE SETS LARGELY OVERLAP, THE 7.337x IS A MEMORISATION NUMBER, NOT A GENERALISATION ONE,
and section 34's comparison against the editor is unfair -- the editor arm was never scored
by a model that had seen its sites.

Measure the overlap exactly, by (chrom, pos). No inference, just the intersection.
"""
import numpy as np
FEAT="/data/a3a/feat"; PC="/data/a3a/pcawg"
CH=[str(i) for i in range(1,23)]+["X"]
d=np.load(f"{FEAT}/a3a_trainset_v5.npz",allow_pickle=True)
y=d["y"].astype(bool)
tr_chrom=d["chrom"][y]; tr_pos=d["pos"][y].astype(np.int64)
train_donors=set(np.unique(d["donor"]).tolist())
print(f"  v5 positives: {len(tr_pos):,}   training donors: {len(train_donors)}")

s=np.load(f"{PC}/snvs.npz",allow_pickle=True)
keep=np.isin(s["donor"],list(train_donors))
mu_chrom=s["chrom"][keep]; mu_pos=s["pos"][keep].astype(np.int64)
print(f"  PCAWG mutations from those donors: {keep.sum():,}")

# exact intersection on (chrom,pos), chromosome by chromosome
inter=0; on_univ=0
for c in CH:
    tp=np.sort(np.unique(tr_pos[tr_chrom==c]))
    mp=np.unique(mu_pos[mu_chrom==c])
    U=np.load(f"{FEAT}/universe_chr{c}.npz")
    upos=U["pos1"] if "pos1" in U.files else U["pos"].astype(np.int64)+1
    # mutations that land on the universe (the denominator s14 used)
    i=np.searchsorted(upos,mp); i=np.clip(i,0,len(upos)-1)
    onu=mp[upos[i]==mp]; on_univ+=len(onu)
    # of those, how many are v5 training positives
    j=np.searchsorted(tp,onu); j=np.clip(j,0,max(len(tp)-1,0))
    if len(tp): inter+=int((tp[j]==onu).sum())
    del U
print(f"\n  mutations landing on the TCW universe (s14's denominator): {on_univ:,}")
print(f"  of those, ALSO v5 TRAINING POSITIVES:                       {inter:,}")
print(f"  IN-SAMPLE FRACTION OF THE SCORED POSITIVES:                 {100*inter/max(on_univ,1):.1f}%")
print()
if inter/max(on_univ,1) > 0.5:
    print("  *** SECTION 34 IS SUBSTANTIALLY IN-SAMPLE. The 7.337x is a memorisation number.")
    print("  *** The comparison against the editor arm is UNFAIR and must be redone with a")
    print("  *** model that has not seen the sites it scores.")
else:
    print("  the scored set is mostly NOT training data; section 34 survives on this axis")
