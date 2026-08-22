#!/usr/bin/env python
"""QA checks 1/2/4 on the NEW artifact this tick: a3a_trainset_v5cov.npz.

A trainset I built myself gets the same leak checks as one a pipeline built. The builder
verified trinuc and strand matching as it wrote; that is the author checking his own work,
so it is re-derived here from the file on disk instead.

Sampling is RANDOM, not a head slice -- these files are block-ordered and a head slice
measures write order (that mistake cost a false alarm two ticks ago)."""
import numpy as np
from collections import Counter

FEAT = "/data/a3a/feat"
B = np.array(list("ACGT"))
rng = np.random.default_rng(0)

for name in ("a3a_trainset_v5.npz", "a3a_trainset_v5cov.npz"):
    d = np.load(f"{FEAT}/{name}", allow_pickle=True)
    y, tri, strand, win = d["y"], d["tri"], d["strand"], d["win"]
    mid = win.shape[1] // 2
    N = len(y)
    i = rng.choice(N, min(400000, N), replace=False); i.sort()
    Y, T, S, W = y[i], tri[i], strand[i], win[i]
    P, Ng = Y == 1, Y == 0
    br = (y == 1).mean()
    print(f"\n=== {name}  n={N:,} pos={int((y==1).sum()):,} "
          f"base rate={br:.6f}  ARITHMETIC CEILING={1/br:.3f}x ===")

    # CHECK 1a focal base
    fp = Counter(B[W[P, mid]].tolist()); fn = Counter(B[W[Ng, mid]].tolist())
    print(f"  focal C: pos {fp['C']/sum(fp.values()):.6f}  neg {fn['C']/sum(fn.values()):.6f}")

    # CHECK 1b flanks at -1/+1 must be identical between classes
    for off in (-1, 1):
        a = Counter(B[W[P, mid+off]].tolist()); b = Counter(B[W[Ng, mid+off]].tolist())
        ta, tb = sum(a.values()), sum(b.values())
        ks = sorted(set(a) | set(b))
        dm = max(abs(a.get(k,0)/ta - b.get(k,0)/tb) for k in ks)
        print(f"  offset {off:+d}: " + " ".join(f"{k} {a.get(k,0)/ta:.4f}/{b.get(k,0)/tb:.4f}" for k in ks)
              + f"  maxdiff={dm:.5f} {'OK' if dm < 0.005 else '*** MISMATCH ***'}")

    # CHECK 1c/1d strand + trinuc
    sv = sorted(set(np.unique(S).tolist()))
    print("  strand: " + "  ".join(f"{v}: {float((S[P]==v).mean()):.4f}/{float((S[Ng]==v).mean()):.4f}" for v in sv))
    tp, tn = Counter(T[P].tolist()), Counter(T[Ng].tolist())
    ta, tb = sum(tp.values()), sum(tn.values())
    dm = max(abs(tp.get(k,0)/ta - tn.get(k,0)/tb) for k in set(tp)|set(tn))
    print("  trinuc: " + "  ".join(f"{k} {tp.get(k,0)/ta:.4f}/{tn.get(k,0)/tb:.4f}" for k in sorted(tp))
          + f"  maxdiff={dm:.5f} {'OK' if dm < 0.005 else '*** MISMATCH ***'}")

    # CHECK 2 YTCA derived
    tcw = np.isin(T, ["TCA", "TCT"])
    for lab, m in (("pos", P & tcw), ("neg", Ng & tcw)):
        yv = B[W[m, mid-2]]
        print(f"  YTCA at TCW [{lab}]: {np.isin(yv, ['C','T']).mean():.4f}  n={int(m.sum()):,}"
              + ("   (locked background 0.6057)" if lab == "neg" else ""))

    # CHECK 4 GC deciles -- positive rate must not slide with local GC
    gc = d["local_gc"][i]
    qs = np.unique(np.percentile(gc, np.arange(0, 101, 10)))
    print("  GC-decile positive-rate ratio: ", end="")
    out = []
    for k in range(len(qs)-1):
        m = (gc >= qs[k]) & (gc < qs[k+1]) if k < len(qs)-2 else (gc >= qs[k]) & (gc <= qs[k+1])
        if m.sum() < 500: continue
        out.append(P[m].mean() / P.mean())
    print(" ".join(f"{v:.3f}" for v in out))
    print(f"    spread {min(out):.3f}-{max(out):.3f} (v3's coverage gradient was 0.621-1.322 for scale)")
