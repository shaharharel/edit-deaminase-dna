#!/usr/bin/env python
"""Are the six hairpin features computed IDENTICALLY for positives and negatives?

This is the check the whole architecture result now rests on. Two independent learners agree
at ~4.9-5.0x on those six features while the best hand rule reaches 2.988x -- a large gap. If
a feature were computed even slightly differently for the two classes, BOTH learners would
inherit the same leak and their agreement would prove nothing. That is precisely bug 2's
shape: a -2 base taken through a ref=='C' branch, correct for one class and not the other.

Method, convention-free on purpose: I do not have the original scanner's exact stem/loop
convention, so instead of reimplementing it (and possibly reproducing its bug) I scan each
window for the longest inverted repeat around the focal base with completely separate code,
then ask whether the relationship between MY number and the STORED number is the same in both
classes. A leak shows up as a class-dependent relationship, whatever the convention is.
"""
import numpy as np

FEAT = "/data/a3a/feat"
d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
W = np.load(f"{FEAT}/a3a_trainset_v5_win1k.npz", allow_pickle=True)
assert (W["pos"] == d["pos"]).all(), "win1k row order differs from v5"
y, stem, loop = d["y"], d["stem"], d["loop"]
win = W["win1k"]; MID = win.shape[1] // 2
N = len(y)

rng = np.random.default_rng(0)
SUB = rng.choice(N, 60000, replace=False); SUB.sort()
print(f"random sample n={len(SUB):,}  base rate {y[SUB].mean():.5f} "
      f"({'OK' if abs(y[SUB].mean()-1/11) < 0.005 else '*** NOT REPRESENTATIVE ***'})")

COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)   # A<->T, C<->G, N->N

def longest_ir(seq, centre, max_loop=12, max_stem=15):
    """Longest inverted repeat whose loop contains the centre. Deliberately naive."""
    best = 0
    for l in range(3, max_loop + 1):
        for off in range(0, l):
            lo = centre - off - 1          # last base before the loop
            hi = centre - off + l          # first base after the loop
            s = 0
            while s < max_stem and lo - s >= 0 and hi + s < len(seq):
                a, b = seq[lo - s], seq[hi + s]
                if a > 3 or b > 3 or COMP[a] != b:
                    break
                s += 1
            if s > best:
                best = s
    return best

mine = np.zeros(len(SUB), np.int16)
for i, r in enumerate(SUB):
    mine[i] = longest_ir(win[r], MID)
    if i % 20000 == 0 and i:
        print(f"  scanned {i:,}", flush=True)

ys = y[SUB]; st = stem[SUB]
print(f"\n=== my independent scan vs the stored 'stem', BY CLASS ===")
print(f"  {'class':>10} {'n':>8} {'mean stored':>12} {'mean mine':>10} {'corr':>7} "
      f"{'mean(mine-stored)':>18}")
for lab, m in (("positives", ys == 1), ("negatives", ys == 0)):
    a, b = st[m].astype(float), mine[m].astype(float)
    c = float(np.corrcoef(a, b)[0, 1])
    print(f"  {lab:>10} {int(m.sum()):>8,} {a.mean():>12.4f} {b.mean():>10.4f} {c:>7.4f} "
          f"{(b - a).mean():>18.4f}")
dc = abs(float(np.corrcoef(st[ys == 1].astype(float), mine[ys == 1])[0, 1]) -
         float(np.corrcoef(st[ys == 0].astype(float), mine[ys == 0])[0, 1]))
dd = abs((mine[ys == 1] - st[ys == 1]).mean() - (mine[ys == 0] - st[ys == 0]).mean())
print(f"\n  class difference in correlation: {dc:.4f}")
print(f"  class difference in mean offset:  {dd:.4f}")
print(f"  {'CONSISTENT -- no class-dependent computation' if dc < 0.02 and dd < 0.05 else '*** CLASS-DEPENDENT: investigate ***'}")

print(f"\n=== and does MY feature reproduce the enrichment? (a leak would not transfer) ===")
base = ys.mean()
for name, s in (("stored stem", st.astype(float)), ("my independent scan", mine.astype(float))):
    sc = s + rng.random(len(s)) * 1e-6
    k = max(int(round(len(s) * 0.001)), 1)
    idx = np.argpartition(-sc, k - 1)[:k]
    print(f"  {name:22s} top0.1% {ys[idx].mean()/base:6.3f}x  (n_top={k})")
