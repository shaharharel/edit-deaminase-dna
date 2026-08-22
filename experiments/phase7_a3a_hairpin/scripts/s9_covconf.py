#!/usr/bin/env python
"""6(d): COVERAGE / MAPPABILITY — the last uncontrolled Phase-1-6 confound.

Negatives are matched on donor + trinucleotide + strand, but NOT on how well short reads
map to the site. PCAWG positives can only be CALLED where the donor's WGS had usable depth
and unique mapping. If hairpin stems are systematically better- or worse-mapped than flat
sequence, the entire 4.96x tail could be a mappability gradient wearing a hairpin costume.

INSTRUMENT: we do not have PCAWG per-donor depth, but we have something nearly as good and
entirely independent -- our own HEK293T Parent WGS, aligned to the same hg19, with measured
per-site coverage at exactly these coordinates (feat/counts_Parent_chr*.npz). Depth there is
driven by mappability and copy number, which is the confounding axis. It is a PROXY: it
carries hg19 mappability faithfully and PCAWG's donor-specific depth not at all.

ORDER OF WORK (deliberate): first ask whether the confound EXISTS -- does coverage differ
between positives and negatives, and does it correlate with stem? Only if it does is there
anything to re-match and retrain for. Reporting "we controlled for X" when X was never
imbalanced is a way of sounding rigorous while learning nothing.
"""
import numpy as np
from collections import defaultdict

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]

d = np.load(f"{FEAT}/a3a_trainset_v3.npz", allow_pickle=True)
chrom, pos, y, stem, tri = d["chrom"], d["pos"], d["y"], d["stem"], d["tri"]
print(f"v3: n={len(y):,}  pos={int((y==1).sum()):,}")

cov = np.full(len(y), -1, dtype=np.int32)
hit = 0
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    # BUG 3 lives here: universe "pos" is 0-BASED, trainset "pos" is 1-BASED, same name.
    # The universe carries an explicitly-named pos1 for exactly this reason -- use it, and
    # verify rather than trust: a wrong offset gives 0% (caught) or ~12% (would not be).
    upos = U["pos1"] if "pos1" in U.files else U["pos"] + 1
    m = chrom == c
    if not m.any():
        continue
    idx = np.searchsorted(upos, pos[m])
    idx = np.clip(idx, 0, len(upos) - 1)
    ok = upos[idx] == pos[m]
    sub = np.flatnonzero(m)[ok]
    cov[sub] = P["cov"][idx[ok]]
    hit += int(ok.sum())
    del U, P

print(f"joined to Parent coverage: {hit:,} / {len(y):,} ({100*hit/len(y):.2f}%)")
print("  (misses are sites outside the TCW universe build or on unplaced contigs)")
if hit < 0.5 * len(y):
    raise SystemExit(f"*** JOIN FAILED ({100*hit/len(y):.2f}%) -- coordinate convention wrong, "
                     "not a biology result. Refusing to report numbers from a broken join. ***")

j = cov >= 0
Y, C, S = y[j], cov[j], stem[j]
P, N = Y == 1, Y == 0
print(f"\n=== A. does coverage differ between classes at all? ===")
for lab, m in (("positives", P), ("negatives", N)):
    q = np.percentile(C[m], [10, 25, 50, 75, 90])
    print(f"  {lab:10s} n={int(m.sum()):>10,}  mean={C[m].mean():7.2f}  p10/25/50/75/90 = "
          + "/".join(f"{x:.0f}" for x in q))
delta = C[P].mean() - C[N].mean()
pooled_sd = np.sqrt((C[P].var() + C[N].var()) / 2)
print(f"  mean difference {delta:+.2f} reads  =  {delta/pooled_sd:+.4f} pooled SD  "
      + ("(negligible)" if abs(delta/pooled_sd) < 0.05 else "*** IMBALANCED ***"))

print(f"\n=== B. is coverage correlated with the hairpin feature? (the thing that would make it a confound) ===")
for s in (0, 4, 6, 8):
    m = S >= s
    print(f"  stem>={s}: n={int(m.sum()):>10,}  mean cov={C[m].mean():7.2f}  "
          f"frac cov<10 = {np.mean(C[m] < 10):.4f}")
r = np.corrcoef(S.astype(float), C.astype(float))[0, 1]
print(f"  Pearson r(stem, coverage) = {r:+.4f}")

print(f"\n=== C. positive rate BY COVERAGE DECILE (flat = no confound; the GC-decile test's twin) ===")
qs = np.percentile(C, np.arange(0, 101, 10))
base = P.mean()
print(f"  overall positive rate {base:.6f} (1/11 = {1/11:.6f})")
for i in range(10):
    lo, hi = qs[i], qs[i + 1]
    m = (C >= lo) & (C <= hi) if i == 9 else (C >= lo) & (C < hi)
    if m.sum() < 1000:
        continue
    pr = P[m].mean()
    print(f"  decile {i:2d}  cov {lo:5.0f}-{hi:5.0f}  n={int(m.sum()):>9,}  pos rate {pr:.6f}  "
          f"ratio {pr/base:6.3f}x  mean stem {S[m].mean():.3f}")
