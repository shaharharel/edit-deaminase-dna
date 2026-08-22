#!/usr/bin/env python
"""6(d) part 3: build a COVERAGE-MATCHED version of v5.

Part 2 showed the hairpin axis survives coverage stratification (MH 1.323 vs crude 1.349 at
stem>=8). That clears the feature, not the trained model: the model reads 81 bp of sequence
and could exploit repeat/mappability-correlated motifs that stratification on one feature
does not touch. The only way to test the MODEL is to remove the imbalance from its training
data and retrain.

Matching is on (trinucleotide, strand, coverage decile). Trinuc and strand were already
matched -- adding coverage to the same key preserves that, it does not replace it.

Coverage comes from our Parent HEK293T WGS at the same hg19 coordinates. The universe's
"pos" is 0-BASED and the trainset's "pos" is 1-BASED under the same field name (bug 3), so
this uses the explicit pos1 field and aborts if the join is incomplete.

Negatives are subsampled without replacement; positives are untouched. The achieved base
rate is reported, because the arithmetic ceiling for the enrichment metric is 1/base_rate
and it must be recomputed if matching moves it off 1/11.
"""
import numpy as np
from collections import defaultdict

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
SRC, DST = "a3a_trainset_v5.npz", "a3a_trainset_v5cov.npz"
rng = np.random.default_rng(0)

d = np.load(f"{FEAT}/{SRC}", allow_pickle=True)
chrom, pos, y, tri, strand = d["chrom"], d["pos"], d["y"], d["tri"], d["strand"]
print(f"{SRC}: n={len(y):,} pos={int((y==1).sum()):,} base rate={(y==1).mean():.6f}")

cov = np.full(len(y), -1, dtype=np.int32)
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    upos = U["pos1"] if "pos1" in U.files else U["pos"] + 1
    m = chrom == c
    if not m.any(): continue
    idx = np.clip(np.searchsorted(upos, pos[m]), 0, len(upos) - 1)
    ok = upos[idx] == pos[m]
    cov[np.flatnonzero(m)[ok]] = P["cov"][idx[ok]]
    del U, P
frac = float((cov >= 0).mean())
print(f"joined to Parent coverage: {frac:.4f}")
if frac < 0.99:
    raise SystemExit(f"*** JOIN INCOMPLETE ({frac:.4f}) -- coordinate convention, not biology ***")

P_ = y == 1
edges = np.unique(np.percentile(cov[P_], np.arange(0, 101, 10)))
cbin = np.clip(np.searchsorted(edges, cov, side="right") - 1, 0, len(edges) - 2)
print(f"coverage bins from the POSITIVE distribution: {[int(e) for e in edges]}")

print("\nBEFORE matching:")
print(f"  mean cov  pos {cov[P_].mean():7.2f}   neg {cov[~P_].mean():7.2f}   "
      f"delta {cov[P_].mean()-cov[~P_].mean():+.2f}")

key_p = defaultdict(int)
idx_n = defaultdict(list)
for i in range(len(y)):
    k = (tri[i], int(strand[i]), int(cbin[i]))
    if y[i] == 1: key_p[k] += 1
    else: idx_n[k].append(i)

keep = [np.flatnonzero(P_)]
short = tot_want = tot_got = 0
for k, npos in key_p.items():
    want = npos * 10
    avail = idx_n.get(k, [])
    tot_want += want
    if len(avail) <= want:
        take = np.array(avail, dtype=np.int64); short += 1
    else:
        take = rng.choice(np.array(avail, dtype=np.int64), want, replace=False)
    tot_got += len(take)
    keep.append(take)
sel = np.sort(np.concatenate(keep))
print(f"  strata: {len(key_p):,}   negative-short strata: {short:,}   "
      f"negatives wanted {tot_want:,} got {tot_got:,} ({tot_got/max(tot_want,1):.4f})")

ys = y[sel]; cs = cov[sel]; ps = ys == 1
print("\nAFTER matching:")
print(f"  n={len(sel):,}  pos={int(ps.sum()):,}  base rate={ps.mean():.6f}  "
      f"(ceiling {1/ps.mean():.3f}x)")
print(f"  mean cov  pos {cs[ps].mean():7.2f}   neg {cs[~ps].mean():7.2f}   "
      f"delta {cs[ps].mean()-cs[~ps].mean():+.2f}")
for q in (10, 50, 90):
    print(f"  p{q:<2d} cov  pos {np.percentile(cs[ps], q):7.1f}   neg {np.percentile(cs[~ps], q):7.1f}")
# trinuc+strand must still be matched -- adding a key cannot be allowed to break the old one
for name, arr in (("trinuc", tri[sel]), ("strand", strand[sel].astype(str))):
    a = {k: (arr[ps] == k).mean() for k in np.unique(arr)}
    b = {k: (arr[~ps] == k).mean() for k in np.unique(arr)}
    dm = max(abs(a[k] - b.get(k, 0)) for k in a)
    print(f"  {name} still matched: maxdiff {dm:.5f} " + ("OK" if dm < 0.01 else "*** BROKEN ***"))

out = {k: d[k][sel] for k in d.files}
out["cov"] = cs
np.savez_compressed(f"{FEAT}/{DST}", **out)
print(f"\nwrote {DST}")
