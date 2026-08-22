#!/usr/bin/env python
"""Does the tail need a MODEL at all, or is it a threshold?

struct_only -- six hairpin-geometry features through an MLP -- reached 4.881x, and the
foundation models reached 1.1x. If a hand-written rule on ONE of those features gets the same
4.9x, then the deliverable for a regulatory-grade safety gate is a transparent, auditable
rule, not a network. A reviewer can check "stem length >= k"; nobody can check 100M weights.

No training and no folds are needed: a fixed rule has nothing to fit, so it cannot overfit,
and its enrichment is directly comparable to a cross-validated model's pooled-OOF number
(if anything the comparison FAVOURS the model, which got to see 80% of the labels).

Random baseline beside every number. Base rate 1/11, so 11.000x is the arithmetic ceiling.
"""
import numpy as np

FEAT = "/data/a3a/feat"
d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
y = d["y"].astype(np.float64); N = len(y); base = y.mean()
print(f"v5 n={N:,} base rate {base:.6f} ceiling {1/base:.3f}x\n")

rng = np.random.default_rng(0)
def tail(score, pct=0.1, jitter=True):
    s = score.astype(np.float64)
    if jitter:                       # break ties at random, not by row order
        s = s + rng.random(N) * 1e-6
    k = max(int(round(N * pct / 100)), 1)
    idx = np.argpartition(-s, k - 1)[:k]
    return y[idx].mean() / base, k

RULES = {
    "stem length":                d["stem"],
    "hp_score":                   d["hp_score"],
    "stem, tie-break hp_score":   d["stem"].astype(np.float64) + 1e-3 * d["hp_score"],
    "gc_pairs":                   d["gc_pairs"],
    "loop length (negated)":      -d["loop"].astype(np.float64),
    "local GC":                   d["local_gc"],
}
print(f"  {'rule (no training, no folds)':32s} {'top0.1%':>9} {'top1%':>8} {'random':>8}")
for name, s in RULES.items():
    e01, k = tail(np.asarray(s), 0.1)
    e1, _ = tail(np.asarray(s), 1.0)
    r, _ = tail(rng.permutation(N).astype(np.float64), 0.1, jitter=False)
    print(f"  {name:32s} {e01:8.3f}x {e1:7.3f}x {r:7.3f}x")

print(f"\n  learned, for comparison (pooled OOF, 5 folds by held-out chromosome, n_top=924):")
print(f"    {'GB one-hot+hairpin':32s} {5.280:8.3f}x")
print(f"    {'MLP struct_only (6 feat)':32s} {4.881:8.3f}x")
print(f"    {'MLP ntv2_only (1024 feat)':32s} {1.155:8.3f}x")

# how many sites even reach the top stem values? the rule's resolution is bounded by ties
st = d["stem"]
print(f"\n  stem-length distribution (why the rule saturates):")
for k in sorted(set(st.tolist()))[-8:]:
    n = int((st >= k).sum())
    print(f"    stem >= {k:2d}: {n:>9,} sites ({100*n/N:6.3f}%)  positive rate ratio "
          f"{(y[st >= k].mean()/base if n else float('nan')):.3f}x")
print(f"  top 0.1% needs {int(round(N*0.001)):,} sites, so any threshold rule must cut inside")
print(f"  a tied block -- which is why the tie-break matters and why a model can beat it.")
