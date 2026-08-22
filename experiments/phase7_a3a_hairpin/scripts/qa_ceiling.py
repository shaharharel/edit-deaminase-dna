#!/usr/bin/env python
"""Is the ARCHITECTURE or the DATA the binding constraint on the 5.280x tail?

Before spending a night on encoders it is worth asking what the attainable maximum even is.
The arithmetic ceiling is 11.000x (base rate 1/11), but that assumes a perfect ranker AND
that the labels mark every site where the process can act. Neither holds:

  - a site is labelled positive because SOME donor happened to mutate it. A site that is
    equally susceptible but was not hit in these 133 donors is labelled NEGATIVE.
  - so the negative set is contaminated with would-be positives, and no model, of any
    architecture, can rank them correctly. That caps the achievable tail well below 11x.

The measurable proxy is DONOR SPLIT-HALF REPRODUCIBILITY. Split the donors into two disjoint
halves. Sites positive in half A are, for half B, exactly "susceptible but not observed".
The overlap between the halves -- against what chance predicts -- estimates how much of the
signal is site-intrinsic and therefore learnable at all.

If sites are essentially never shared between disjoint donor halves, then site identity is
not predictable even in principle, and a model can only predict a PROPENSITY. The ceiling is
then set by how much propensity varies across sites -- which is what the hairpin feature
already captures, and would mean 5.280x is close to the data's limit rather than the
architecture's.
"""
import numpy as np
from collections import Counter

FEAT = "/data/a3a/feat"
d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
y, chrom, pos, don, stem = d["y"], d["chrom"], d["pos"], d["donor"], d["stem"]
P = y == 1
key = np.char.add(np.char.add(chrom.astype(str), ":"), pos.astype(str))
kp, dp = key[P], don[P].astype(str)
print(f"v5 positives: {len(kp):,}  distinct sites: {len(set(kp)):,}  donors: {len(set(dp)):,}")
print("NOTE: v5 holds a 84k SAMPLE of the 2.38M available positives, so absolute overlap")
print("here understates the full cohort. The obs/exp ratio is the informative quantity.")

c = Counter(kp.tolist())
mult = Counter(c.values())
print(f"\n=== donor recurrence of a site (within this sample) ===")
for k in sorted(mult)[:5]:
    print(f"  seen in {k} donor-record(s): {mult[k]:,} sites ({100*mult[k]/len(c):.3f}%)")

donors = sorted(set(dp.tolist()))
rng = np.random.default_rng(0)
perm = rng.permutation(len(donors))
A = {donors[i] for i in perm[:len(donors)//2]}
inA = np.array([x in A for x in dp])
sA, sB = set(kp[inA].tolist()), set(kp[~inA].tolist())
ov = sA & sB
# chance expectation over the eligible TCW universe the negatives were drawn from
UNIV = len(set(key.tolist())) / max((y == 0).mean(), 1e-9)   # crude scale of the site space
exp = len(sA) * len(sB) / UNIV
print(f"\n=== disjoint donor halves ({len(A)} vs {len(donors)-len(A)} donors) ===")
print(f"  sites in half A: {len(sA):,}   half B: {len(sB):,}")
print(f"  shared: {len(ov):,}   expected by chance over ~{UNIV:,.0f} sites: {exp:,.1f}")
print(f"  obs/exp: {len(ov)/max(exp,1e-9):.2f}x")
print(f"  reproducible fraction of half A: {len(ov)/max(len(sA),1):.5f}")

print(f"\n=== is the SHARED fraction more hairpin-rich than the private one? ===")
hp = {k: [] for k in ("shared", "private")}
for k, s in zip(kp.tolist(), stem[P].tolist()):
    hp["shared" if k in ov else "private"].append(s)
for lab in ("shared", "private"):
    a = np.array(hp[lab])
    if len(a) < 20:
        print(f"  {lab:8s} n={len(a):,} (too few)"); continue
    print(f"  {lab:8s} n={len(a):>7,}  stem>=6 {np.mean(a>=6):.4f}  stem>=8 {np.mean(a>=8):.4f}")
bg = stem[~P]
print(f"  {'negatives':8s} n={len(bg):>7,}  stem>=6 {np.mean(bg>=6):.4f}  stem>=8 {np.mean(bg>=8):.4f}")
print("\n  If shared sites are far more hairpin-rich, the learnable part of the signal is the")
print("  hairpin part -- and a better sequence encoder has little left to find.")
