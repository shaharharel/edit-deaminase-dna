#!/usr/bin/env python
"""QA: are clone2 (1.069) and clone5 (1.323) at stem 8 actually DIFFERENT?

Two point estimates read against two SEPARATE nulls is not a comparison. n is 43 and
53 - the regime where this project's A3A-vs-A3B claim died at n=97. Before treating
"stem 8 is unresolved between clones" as a finding, test whether the two clones are
statistically distinguishable at all, and what the POOLED estimate looks like.

Recomputes the s7c endpoint-B components directly from counts_:
    silent = Parent alt==0 & Parent cov>=15          (the germline mask s7c uses)
    elig   = silent & sample cov>=8
    spec   = elig & sample alt>=2
    enr    = P(hairpin | spec) / P(hairpin | elig)
Then: a 2x2 test between clones on (spec & hairpin) vs (spec & not hairpin), and a
pooled enrichment with its own binomial null.
"""
import numpy as np
from scipy import stats

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
STEMS = [6, 7, 8]
EDITORS = ["P66-A3A-Y130F-clone2", "P66-A3A-Y130F-clone5"]
CAL = "P66-D10A-clone1"

acc = {s: {k: np.zeros(4, dtype=np.int64) for k in EDITORS + [CAL]} for s in STEMS}
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    stem = U["stem"]
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    silent = (P["alt"] == 0) & (P["cov"] >= 15)
    for name in EDITORS + [CAL]:
        D = np.load(f"{FEAT}/counts_{name}_chr{c}.npz")
        elig = silent & (D["cov"] >= 8)
        spec = elig & (D["alt"] >= 2)
        for s in STEMS:
            h = stem >= s
            a = acc[s][name]
            a[0] += int((spec & h).sum())   # spec & hairpin
            a[1] += int(spec.sum())         # spec total
            a[2] += int((elig & h).sum())   # elig & hairpin
            a[3] += int(elig.sum())         # elig total

rng = np.random.default_rng(20260821)
print("=== endpoint-B components recomputed from counts_ (independent of s7c) ===")
for s in STEMS:
    print(f"\n--- stem >= {s}")
    print(f"  {'sample':24s} {'sh':>6s} {'sn':>9s} {'p_bg':>8s} {'enr':>7s}")
    vals = {}
    for name in EDITORS + [CAL]:
        sh, sn, bh, bn = acc[s][name]
        pb = bh / max(bn, 1)
        enr = (sh / max(sn, 1)) / max(pb, 1e-12)
        vals[name] = (sh, sn, pb, enr)
        print(f"  {name:24s} {sh:>6,} {sn:>9,} {pb:>8.5f} {enr:>6.3f}")
    # between-clone 2x2 on the two editors
    (s2, n2, _, e2) = vals[EDITORS[0]]
    (s5, n5, _, e5) = vals[EDITORS[1]]
    tab = [[s2, n2 - s2], [s5, n5 - s5]]
    odds, p = stats.fisher_exact(tab)
    print(f"  BETWEEN CLONES  clone2 enr={e2:.3f} (n_hp={s2})  clone5 enr={e5:.3f} (n_hp={s5})")
    print(f"                  Fisher exact on hairpin share: p = {p:.4f}  "
          f"{'-> NOT distinguishable' if p > 0.05 else '-> distinguishable'}")
    # pooled editor estimate
    sh_p, sn_p = s2 + s5, n2 + n5
    pb_p = (vals[EDITORS[0]][2] * n2 + vals[EDITORS[1]][2] * n5) / (n2 + n5)
    enr_p = (sh_p / sn_p) / pb_p
    draws = rng.binomial(sn_p, pb_p, 4000) / sn_p / pb_p
    lo, hi = np.percentile(draws, [2.5, 97.5])
    pv = float((draws >= enr_p).mean())
    print(f"  POOLED EDITORS  enr={enr_p:.3f}  null95 {lo:.2f}-{hi:.2f}  p={pv:.4f}  "
          f"n_hp={sh_p}  (calibrator {vals[CAL][3]:.3f})")
