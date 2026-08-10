#!/usr/bin/env python
"""
S7b — editor-specific hairpin enrichment with a proper null DISTRIBUTION.

s7_editor_test.py reports a single binomial draw as its RANDOM baseline. That
shows the null is "about 1.0" but not whether an observed 1.53x lies outside it.
This reports:

  * the null 95% interval from 2000 draws, and an empirical one-sided p
  * a CI on the OBSERVED enrichment -- the high-stem bins carry only a few
    hundred sites, and that uncertainty is usually the binding constraint

Both are needed for the actual question: the pre-registered bar is a >=3x
increment over a deaminase-free baseline of 1.28-1.63x, so an editor must clear
~4.6-4.9x. Point estimates on either side cannot settle that.

Same mandatory germline filter as s7 (see HANDOFF.md discipline #6).

Usage: s7b_stats.py <editor> <control1,control2,...|''> [parent]
"""
import sys, glob, re
import numpy as np

FEAT = "/mnt/data/a3a/feat"
editor = sys.argv[1]
controls = [c for c in sys.argv[2].split(",") if c] if len(sys.argv) > 2 else []
parent = sys.argv[3] if len(sys.argv) > 3 else "Parent"
rng = np.random.default_rng(20260810)
STEMS = [4, 5, 6, 7, 8]
NDRAW = 2000


def chroms(s):
    return {re.search(r"_chr([0-9XY]+)\.npz", f).group(1)
            for f in glob.glob(f"{FEAT}/counts_{s}_chr*.npz")}


common = chroms(editor)
for c in controls + [parent]:
    common &= chroms(c)
common = sorted(common, key=lambda c: (len(c), c))
print(f"editor={editor}  controls={controls}  parent={parent}  chroms={len(common)}")

n_spec = n_bg = 0
hp_spec = {s: 0 for s in STEMS}
hp_bg = {s: 0 for s in STEMS}
for c in common:
    stem = np.load(f"{FEAT}/universe_chr{c}.npz")["stem"]
    E = np.load(f"{FEAT}/counts_{editor}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{parent}_chr{c}.npz")
    keep = (P["alt"] == 0) & (P["cov"] >= 15)
    for cs in controls:
        C = np.load(f"{FEAT}/counts_{cs}_chr{c}.npz")
        keep &= (C["alt"] == 0) & (C["cov"] >= 8)
    keep &= E["cov"] >= 8
    spec = keep & (E["alt"] >= 2)
    n_spec += int(spec.sum()); n_bg += int(keep.sum())
    for s in STEMS:
        h = stem >= s
        hp_spec[s] += int((spec & h).sum())
        hp_bg[s] += int((keep & h).sum())

print(f"eligible background: {n_bg:,}    editor-specific: {n_spec:,}\n")
print(f"{'stem':>5} {'enrich':>8} {'obs 95%CI':>14} {'null 95%CI':>14} {'p':>8} {'n_hp':>7}")
for s in STEMS:
    k = hp_spec[s]
    ps = k / max(n_spec, 1)
    pb = hp_bg[s] / max(n_bg, 1)
    enr = ps / max(pb, 1e-12)
    draws = rng.binomial(n_spec, pb, NDRAW) / max(n_spec, 1) / max(pb, 1e-12)
    lo_n, hi_n = np.percentile(draws, [2.5, 97.5])
    pval = float((draws >= enr).mean())
    se = 1.0 / np.sqrt(max(k, 1))           # Poisson relative SE on the count
    lo_o, hi_o = enr * (1 - 1.96 * se), enr * (1 + 1.96 * se)
    print(f"{s:>5} {enr:>8.3f} {f'{lo_o:.2f}-{hi_o:.2f}':>14} "
          f"{f'{lo_n:.2f}-{hi_n:.2f}':>14} {pval:>8.4f} {k:>7,}")

print("\nInterpreting: an enrichment is only interesting if its observed CI lies")
print("above the null CI. For an EDITOR claim the bar is higher still -- the")
print("observed CI must clear ~4.6-4.9x, i.e. >=3x over the replicated")
print("deaminase-free baseline of 1.28-1.63x. Not above 1.0.")
