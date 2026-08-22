#!/usr/bin/env python
"""CROSS-FAMILY CONTROL, preview on the one Lj-BE clone whose counts are complete.

Lj-BE carries lamprey CDA1 (PmCDA1), a cytidine deaminase from a different family than
APOBEC3. PRE-REGISTERED PREDICTION, written before looking: if the hairpin preference is an
APOBEC-specific property of the deaminase, Lj-BE should sit near the deaminase-free
calibrator, NOT above it. If Lj-BE shows the same hairpin enrichment as the APOBEC arms,
that argues the enrichment is a property of the DNA substrate (exposed ssDNA in stems) and
not of the enzyme -- which would be a more interesting result and a weaker editor claim.

THIS IS 1 OF 3 CLONES. clone5 and clone12 land ~01:30 and ~02:15 UTC. Numbers here are a
preview with n stated; the pooled 3-clone value is the one that counts. The A3A-Y130F arm
is the cautionary case -- its apparent stem-8 signal came entirely from sites shared between
two clones, invisible until the second clone existed.

Endpoint B only (hairpin share among editor-specific sites). No burden endpoint: the
deaminase-free calibrator lives on the other node.
"""
import numpy as np

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
STEMS = [6, 7, 8]
SAMPLES = ["P66-Lj-BE-clone3", "P66-eA3A-RL1-clone1", "P66-eA3A-RL1-clone2", "P66-eA3A-RL1-clone5"]
rng = np.random.default_rng(0)

bg_h = {s: 0 for s in STEMS}; bg_n = 0
res = {}
for smp in SAMPLES:
    sh = {s: 0 for s in STEMS}; sn = 0
    for c in CH:
        U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]
        P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
        D = np.load(f"{FEAT}/counts_{smp}_chr{c}.npz")
        silent = (P["alt"] == 0) & (P["cov"] >= 15)
        spec = silent & (D["cov"] >= 8) & (D["alt"] >= 2)
        sn += int(spec.sum())
        for s in STEMS:
            sh[s] += int((spec & (stem >= s)).sum())
        if smp == SAMPLES[0]:
            bg_n += int(silent.sum())
            for s in STEMS:
                bg_h[s] += int((silent & (stem >= s)).sum())
    res[smp] = (sh, sn)

print(f"background (Parent-silent, cov>=15) n={bg_n:,}")
for s in STEMS:
    print(f"  p_bg(stem>={s}) = {bg_h[s]/bg_n:.6f}")

print(f"\n{'sample':22s} {'n_spec':>9s} " + " ".join(f"{'stem'+str(s):>22s}" for s in STEMS))
for smp in SAMPLES:
    sh, sn = res[smp]
    cells = []
    for s in STEMS:
        p = bg_h[s] / bg_n
        enr = (sh[s] / sn) / p if sn else float("nan")
        # binomial null on the OBSERVED set size: what range does chance alone give?
        draws = rng.binomial(sn, p, 2000) / sn / p
        lo, hi = np.percentile(draws, [2.5, 97.5])
        star = "*" if (enr < lo or enr > hi) else " "
        cells.append(f"{enr:6.3f}x[{lo:.2f}-{hi:.2f}]{star}n={sh[s]:<6d}")
    print(f"{smp:22s} {sn:>9,} " + " ".join(cells))
print("\n  * = outside the 95% binomial null for a set of that size.")
print("  Read against the DEAMINASE-FREE CALIBRATOR 1.28-1.63x, never against 1.0.")
