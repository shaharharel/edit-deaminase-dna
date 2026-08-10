#!/usr/bin/env python
"""
S7 — editor-specific hairpin enrichment, with the calibrator built in.

Usage:  s7_editor_test.py <editor_sample> <control1,control2,...> [parent_sample]

The germline filter is not optional and not a flag. QA showed 35.5% of raw alt>=2
calls are germline-like, and because germline carries HIGH alt counts it would
dominate the top of any alt-ranked list -- which is exactly the top-K% endpoint.
So an editor-specific site is defined as:

    alt >= 2 AND cov >= 8   in the editor clone
    alt == 0 AND cov >= 15  in Parent/background
    alt == 0 AND cov >= 8   in EVERY deaminase-free control

That removes germline and shared clonal APOBEC3 background in one step (the same
STEP-3 control that salvaged the earlier Doman site-recurrence analysis).

Every enrichment is printed beside a label-shuffled RANDOM baseline, and n is
printed for every effect -- the A3A-vs-A3B claim died from a small-n illusion
(+0.355 at n=53 -> +0.079 at n=97), so effect sizes without n are not reportable.

VALIDATION MODE: run with a deaminase-free clone as the "editor". Both arms are
then deaminase-free, so every enrichment MUST come back ~1.0. If it does not, the
analysis is broken, not the biology.
"""
import sys, glob, re
import numpy as np

FEAT = "/mnt/data/a3a/feat"
editor = sys.argv[1]
controls = [c for c in sys.argv[2].split(",") if c]
parent = sys.argv[3] if len(sys.argv) > 3 else "Parent"
rng = np.random.default_rng(20260810)

STEMS = [4, 5, 6, 7, 8]


def chroms_available(s):
    return {re.search(r"_chr([0-9XY]+)\.npz", f).group(1)
            for f in glob.glob(f"{FEAT}/counts_{s}_chr*.npz")}


common = chroms_available(editor)
for c in controls + [parent]:
    common &= chroms_available(c)
common = sorted(common, key=lambda c: (len(c), c))
if not common:
    print("no chromosomes shared by all samples"); sys.exit(1)
print(f"editor={editor}  controls={controls}  parent={parent}")
print(f"chromosomes usable: {len(common)}  {common}")

n_univ = n_cov = n_spec = 0
hp_spec = {s: 0 for s in STEMS}
hp_bg = {s: 0 for s in STEMS}
n_bg = 0

for c in common:
    u = np.load(f"{FEAT}/universe_chr{c}.npz")
    stem = u["stem"]
    E = np.load(f"{FEAT}/counts_{editor}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{parent}_chr{c}.npz")
    ok_ctrl = (P["alt"] == 0) & (P["cov"] >= 15)
    for cs in controls:
        C = np.load(f"{FEAT}/counts_{cs}_chr{c}.npz")
        ok_ctrl &= (C["alt"] == 0) & (C["cov"] >= 8)
    ecov = E["cov"] >= 8
    spec = ok_ctrl & ecov & (E["alt"] >= 2)
    bg = ok_ctrl & ecov                      # the eligible background: same filter,
                                             # minus the requirement of being edited
    n_univ += len(stem); n_cov += int(ecov.sum())
    n_spec += int(spec.sum()); n_bg += int(bg.sum())
    for s in STEMS:
        h = stem >= s
        hp_spec[s] += int((spec & h).sum())
        hp_bg[s] += int((bg & h).sum())

print(f"\nuniverse sites                 : {n_univ:,}")
print(f"editor cov>=8                  : {n_cov:,}")
print(f"eligible background (filtered) : {n_bg:,}")
print(f"EDITOR-SPECIFIC sites          : {n_spec:,}   "
      f"({n_spec/max(n_bg,1)*1e6:.1f} per Mb of eligible)")

if n_spec == 0:
    print("no editor-specific sites -- nothing to test"); sys.exit(0)

print(f"\n{'stem':>6} {'spec_frac':>10} {'bg_frac':>10} {'enrich':>8} {'RANDOM':>8} {'n_spec':>8}")
for s in STEMS:
    ps = hp_spec[s] / max(n_spec, 1)
    pb = hp_bg[s] / max(n_bg, 1)
    # RANDOM baseline: draw n_spec sites at random from the eligible background
    # and ask how often they are hairpins -- i.e. the enrichment expected by chance.
    draw = rng.binomial(n_spec, pb) / max(n_spec, 1)
    rnd = draw / max(pb, 1e-12)
    print(f"{s:>6} {ps:>10.6f} {pb:>10.6f} {ps/max(pb,1e-12):>8.3f} {rnd:>8.3f} {hp_spec[s]:>8,}")

print("\nn is printed for every effect. Treat any bin with n_spec < ~100 as noise.")
print("If this was run with a deaminase-free clone as 'editor', every enrichment")
print("must be ~1.0 -- anything else is an analysis bug, not biology.")
