#!/usr/bin/env python
"""Before the Lj-BE preview is believed: does clone3's site set look like real clonal
mutation, or like the artefact that disqualified D10A-clone1?

Lj-BE-clone3 carries 1,301,573 Parent-specific sites -- 2.6x eA3A-clone5, 7x eA3A-clone2,
and ~65x the A3A-Y130F clones. That spread is the same shape as D10A-clone1, which had 7.5x
everyone else and turned out to be a sequencing defect: 94.8% of its sites sat at VAF<0.05
and only 1.4% at VAF>=0.35. A clone with a huge low-VAF tail is calling noise, and noise has
its own hairpin behaviour.

Run the same discriminator: the VAF profile, plus depth, plus what the enrichment looks like
when the low-VAF tail is removed. If the hairpin signal only exists in the VAF<0.05 fraction,
the preview means nothing."""
import numpy as np

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
SAMPLES = ["P66-Lj-BE-clone3", "P66-eA3A-RL1-clone1", "P66-eA3A-RL1-clone5"]
BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.35), (0.35, 1.01)]
STEMS = [6, 8]

bg = {s: [0, 0] for s in STEMS}
prof = {s: np.zeros(len(BINS), dtype=np.int64) for s in SAMPLES}
hi_h = {s: {t: 0 for t in STEMS} for s in SAMPLES}; hi_n = {s: 0 for s in SAMPLES}
lo_h = {s: {t: 0 for t in STEMS} for s in SAMPLES}; lo_n = {s: 0 for s in SAMPLES}
depth = {s: [0, 0] for s in SAMPLES}

for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    silent = (P["alt"] == 0) & (P["cov"] >= 15)
    for t in STEMS:
        bg[t][0] += int((silent & (stem >= t)).sum()); bg[t][1] += int(silent.sum())
    for smp in SAMPLES:
        D = np.load(f"{FEAT}/counts_{smp}_chr{c}.npz")
        cov, alt = D["cov"], D["alt"]
        spec = silent & (cov >= 8) & (alt >= 2)
        depth[smp][0] += int(cov.sum()); depth[smp][1] += len(cov)
        v = np.zeros(len(cov)); nz = cov > 0
        v[nz] = alt[nz] / cov[nz]
        vv = v[spec]
        for i, (lo, hi) in enumerate(BINS):
            prof[smp][i] += int(((vv >= lo) & (vv < hi)).sum())
        himask = spec & (v >= 0.15)
        lomask = spec & (v < 0.05)
        hi_n[smp] += int(himask.sum()); lo_n[smp] += int(lomask.sum())
        for t in STEMS:
            hi_h[smp][t] += int((himask & (stem >= t)).sum())
            lo_h[smp][t] += int((lomask & (stem >= t)).sum())

print(f"{'sample':22s} {'mean depth':>10s} {'n_spec':>10s} {'<.05':>7s} {'.05-.15':>8s} {'.15-.35':>8s} {'>=.35':>7s}")
for smp in SAMPLES:
    t = max(prof[smp].sum(), 1)
    print(f"{smp:22s} {depth[smp][0]/depth[smp][1]:>10.2f} {t:>10,} "
          + " ".join(f"{100*x/t:6.1f}%" for x in prof[smp]))
print("  D10A-clone1, DISQUALIFIED for comparison:  94.8% / -- / -- / 1.4%")

print(f"\n=== hairpin enrichment SPLIT BY VAF (does the signal live in the noise tail?) ===")
for t in STEMS:
    p = bg[t][0] / bg[t][1]
    print(f"  stem>={t}  (p_bg={p:.6f})")
    for smp in SAMPLES:
        e_hi = (hi_h[smp][t] / hi_n[smp]) / p if hi_n[smp] else float("nan")
        e_lo = (lo_h[smp][t] / lo_n[smp]) / p if lo_n[smp] else float("nan")
        print(f"    {smp:22s} VAF>=0.15 {e_hi:6.3f}x (n={hi_n[smp]:>8,}, n_hp={hi_h[smp][t]:>6,})   "
              f"VAF<0.05 {e_lo:6.3f}x (n={lo_n[smp]:>9,})")
