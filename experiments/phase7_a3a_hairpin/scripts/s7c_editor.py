#!/usr/bin/env python
"""
S7c — the editor test, both endpoints, with their distinct confounds handled.

Usage: s7c_editor.py <editor> <calibrator> <control1,...|''> [parent]

ENDPOINT A -- BURDEN (editor-specific site count). This is where the >=3x bar
belongs. It is COVERAGE-CONFOUNDED: P(alt>=2 | VAF 0.08) is 0.483 at 20x and
0.841 at 40x, so depth alone inflates burden 1.74x. The P66 editor samples were
sequenced ~1.4x deeper than the PRJNA1042830 calibrators. So burden is reported
WITHIN MATCHED COVERAGE BANDS, never pooled.

ENDPOINT B -- HAIRPIN ENRICHMENT among editor-specific sites. A within-sample
ratio, so largely depth-immune (depth affects numerator and denominator alike).
Quoted against the calibrator's 1.28-1.63x with NO 3x expectation: an A3A-family
editor that targets like endogenous A3A should MATCH the calibrator, not exceed
it. Reported with a 2000-draw null and stratified by trinucleotide context.

The informative pattern for the hypothesis is A high AND B near the calibrator.
B large with A flat is more likely a configuration or filtering artifact.

VALIDATION MODE: pass two deaminase-free clones as editor/calibrator. Endpoint A
must give a burden ratio ~1.0 within each band, and B must match the known
baseline. Anything else is an analysis bug.
"""
import sys, glob, re
import numpy as np

FEAT = "/mnt/data/a3a/feat"
editor, calib = sys.argv[1], sys.argv[2]
controls = [c for c in sys.argv[3].split(",") if c] if len(sys.argv) > 3 else []
parent = sys.argv[4] if len(sys.argv) > 4 else "Parent"
rng = np.random.default_rng(20260810)
BANDS = [(8, 15), (15, 25), (25, 35), (35, 60), (60, 10**6)]
STEMS = [6, 7, 8]


def chroms(s):
    return {re.search(r"_chr([0-9XY]+)\.npz", f).group(1)
            for f in glob.glob(f"{FEAT}/counts_{s}_chr*.npz")}


common = chroms(editor) & chroms(calib) & chroms(parent)
for c in controls:
    common &= chroms(c)
common = sorted(common, key=lambda c: (len(c), c))
if not common:
    print("no shared chromosomes"); sys.exit(1)
print(f"editor={editor}  calibrator={calib}  controls={controls}  parent={parent}")
print(f"chromosomes: {len(common)}\n")

burden = {b: [0, 0, 0, 0] for b in BANDS}      # ed_spec, ed_elig, cal_spec, cal_elig
hp = {w: {s: [0, 0, 0, 0] for s in STEMS} for w in ("ed", "cal")}
ctx_hp = {w: {t: {s: [0, 0, 0, 0] for s in STEMS} for t in (0, 1)} for w in ("ed", "cal")}

for c in common:
    U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]; tri = U["tri"]
    E = np.load(f"{FEAT}/counts_{editor}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{calib}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{parent}_chr{c}.npz")
    silent = (P["alt"] == 0) & (P["cov"] >= 15)
    for cs in controls:
        C = np.load(f"{FEAT}/counts_{cs}_chr{c}.npz")
        silent &= (C["alt"] == 0) & (C["cov"] >= 8)

    for lo, hi in BANDS:
        eb = silent & (E["cov"] >= lo) & (E["cov"] < hi)
        kb = silent & (K["cov"] >= lo) & (K["cov"] < hi)
        b = burden[(lo, hi)]
        b[0] += int((eb & (E["alt"] >= 2)).sum()); b[1] += int(eb.sum())
        b[2] += int((kb & (K["alt"] >= 2)).sum()); b[3] += int(kb.sum())

    for w, D in (("ed", E), ("cal", K)):
        elig = silent & (D["cov"] >= 8)
        spec = elig & (D["alt"] >= 2)
        for s in STEMS:
            h = stem >= s
            a = hp[w][s]
            a[0] += int((spec & h).sum()); a[1] += int(spec.sum())
            a[2] += int((elig & h).sum()); a[3] += int(elig.sum())
            for t in (0, 1):
                m = tri == t
                a2 = ctx_hp[w][t][s]
                a2[0] += int((spec & h & m).sum()); a2[1] += int((spec & m).sum())
                a2[2] += int((elig & h & m).sum()); a2[3] += int((elig & m).sum())

print("ENDPOINT A -- BURDEN, within matched coverage bands (>=3x bar applies here)")
print(f"  {'cov band':>12} {'editor /Mb':>12} {'calib /Mb':>12} {'ratio':>8} {'n_ed':>9}")
for b in BANDS:
    es, ee, ks, ke = burden[b]
    if ee < 1e6 or ke < 1e6:
        continue
    er, kr = es / ee * 1e6, ks / ke * 1e6
    print(f"  {str(b):>12} {er:>12.1f} {kr:>12.1f} {er/max(kr,1e-9):>8.3f} {es:>9,}")

print("\nENDPOINT B -- HAIRPIN ENRICHMENT (quote vs calibrator, NOT vs 1.0)")
print(f"  {'stem':>5} {'editor':>8} {'calib':>8} {'ed null95':>14} {'p_ed':>7} {'n_hp_ed':>8}")
for s in STEMS:
    out = {}
    for w in ("ed", "cal"):
        sh, sn, bh, bn = hp[w][s]
        ps, pb = sh / max(sn, 1), bh / max(bn, 1)
        out[w] = (ps / max(pb, 1e-12), sh, sn, pb)
    e, sh, sn, pb = out["ed"]
    d = rng.binomial(sn, pb, 2000) / max(sn, 1) / max(pb, 1e-12)
    lo, hi = np.percentile(d, [2.5, 97.5]); pv = float((d >= e).mean())
    print(f"  {s:>5} {e:>8.3f} {out['cal'][0]:>8.3f} {f'{lo:.2f}-{hi:.2f}':>14} {pv:>7.4f} {sh:>8,}")

print("\nENDPOINT B by trinucleotide context (mix shift can imitate an editor effect)")
print(f"  {'ctx':>4} {'stem':>5} {'editor':>8} {'calib':>8} {'n_hp_ed':>8}")
for t, lab in [(0, "TCA"), (1, "TCT")]:
    for s in STEMS:
        r = {}
        for w in ("ed", "cal"):
            sh, sn, bh, bn = ctx_hp[w][t][s]
            r[w] = ((sh / max(sn, 1)) / max(bh / max(bn, 1), 1e-12), sh)
        print(f"  {lab:>4} {s:>5} {r['ed'][0]:>8.3f} {r['cal'][0]:>8.3f} {r['ed'][1]:>8,}")

print("\nInformative pattern: A high AND B near the calibrator.")
print("B large with A flat is more likely a configuration/filtering artifact.")
