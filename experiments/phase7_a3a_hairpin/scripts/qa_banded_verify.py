#!/usr/bin/env python
"""Independent re-implementation of the banded endpoint B, to verify the one now used as the
METHOD OF RECORD.

Two reasons this needs checking rather than trusting:
 1. It was written tonight, by me, and it is now the estimator every cross-arm number depends
    on. The patch that inserted it went into a 200-line script through three string anchors.
 2. Its output contains a value that looks wrong: stem>=6 enrichment 0.107 in band (8,15) --
    a tenfold DEPLETION. That is either a real property of low-coverage calling or an
    accumulator bug, and the two are indistinguishable from the number alone.

This recomputes the same quantities from the raw counts with completely separate code, on
chr1 only, and compares. Written to be obviously correct rather than efficient.
"""
import numpy as np

FEAT = "/data/a3a/feat"
ED, CAL, PAR = "P66-Lj-BE-clone3", "P66-eA3A-RL1-clone1", "Parent"
BANDS = [(8, 15), (15, 25), (25, 35), (35, 60), (60, 100), (100, 10**6)]
STEMS = [6, 7, 8]
# chr1 alone showed the same shape as s7c, but "same shape" is not verification. Run all 23
# and compare the MH numbers to s7c's output digit for digit.
CH = [str(i) for i in range(1, 23)] + ["X"]
acc = {b: {s: {"ed": [0, 0, 0, 0], "cal": [0, 0, 0, 0]} for s in STEMS} for b in BANDS}
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem_c = U["stem"]
    P = np.load(f"{FEAT}/counts_{PAR}_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15)
    for lo, hi in BANDS:
        both = sil & (E["cov"] >= lo) & (E["cov"] < hi) & (K["cov"] >= lo) & (K["cov"] < hi)
        if not both.any():
            continue
        for w, D in (("ed", E), ("cal", K)):
            sp = both & (D["cov"] >= 8) & (D["alt"] >= 2)
            for st in STEMS:
                h = stem_c >= st
                a = acc[(lo, hi)][st][w]
                a[0] += int((sp & h).sum()); a[1] += int(sp.sum())
                a[2] += int((both & h).sum()); a[3] += int(both.sum())

print("=== INDEPENDENT re-implementation, all 23 chromosomes ===")
for st in STEMS:
    ne = de = nc = dc = 0.0
    for b in BANDS:
        ed, cal = acc[b][st]["ed"], acc[b][st]["cal"]
        if ed[3] < 500 or ed[1] < 20 or cal[1] < 20:
            continue
        for arr, which in ((ed, "e"), (cal, "c")):
            A = arr[0]; B = arr[1] - arr[0]
            C_ = arr[2] - arr[0]; D_ = arr[3] - arr[1] - C_
            T = A + B + C_ + D_
            if T <= 0: continue
            if which == "e": ne += A * D_ / T; de += B * C_ / T
            else: nc += A * D_ / T; dc += B * C_ / T
    print(f"  stem>={st}  MH editor {ne/de:.3f}   calibrator {nc/dc:.3f}")
print("  s7c reported:  6 -> 1.138 / 1.129    7 -> 1.308 / 1.210    8 -> 1.569 / 1.324")

c = "1"
U = np.load(f"{FEAT}/universe_chr{c}.npz"); stem = U["stem"]
P = np.load(f"{FEAT}/counts_{PAR}_chr{c}.npz")
E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")

silent = (P["alt"] == 0) & (P["cov"] >= 15)
print(f"chr{c}: universe {len(stem):,}  Parent-silent {int(silent.sum()):,}")
print(f"\n  {'band':>14} {'n_both':>10} {'p_bg(s>=6)':>11} {'n_spec_ed':>10} "
      f"{'n_hp_ed':>8} {'enr_ed':>8} {'enr_cal':>8}")
for lo, hi in BANDS:
    both = silent & (E["cov"] >= lo) & (E["cov"] < hi) & (K["cov"] >= lo) & (K["cov"] < hi)
    nb = int(both.sum())
    if nb < 500:
        print(f"  {f'({lo},{hi})':>14} {nb:>10,}  -- too few --"); continue
    h6 = stem >= 6
    pbg = float((both & h6).sum()) / nb
    row = [f"  {f'({lo},{hi})':>14} {nb:>10,} {pbg:>11.5f}"]
    for w, D in (("ed", E), ("cal", K)):
        sp = both & (D["cov"] >= 8) & (D["alt"] >= 2)
        ns = int(sp.sum()); nh = int((sp & h6).sum())
        enr = (nh / ns) / pbg if ns and pbg else float("nan")
        if w == "ed":
            row.append(f" {ns:>10,} {nh:>8,} {enr:>8.3f}")
        else:
            row.append(f" {enr:>8.3f}")
    print("".join(row))

print("\n=== why is the low-coverage band so depleted? VAF is the mechanism, so measure it ===")
for lo, hi in [(8, 15), (35, 60)]:
    both = silent & (E["cov"] >= lo) & (E["cov"] < hi) & (K["cov"] >= lo) & (K["cov"] < hi)
    sp = both & (E["alt"] >= 2)
    cov = E["cov"][sp].astype(float); alt = E["alt"][sp].astype(float)
    v = alt / np.maximum(cov, 1)
    h = stem[sp] >= 6
    if sp.sum() < 50: continue
    print(f"  band ({lo},{hi}): n_called {int(sp.sum()):,}  median VAF {np.median(v):.3f}  "
          f"min VAF possible {2/hi:.3f}")
    if h.sum() > 5:
        print(f"      median VAF at hairpin sites {np.median(v[h]):.3f} (n={int(h.sum())})  "
              f"vs flat {np.median(v[~h]):.3f}")
print("  alt>=2 at cov 8-15 forces VAF>=0.13-0.25; at cov 35-60 it admits VAF>=0.03-0.06.")
print("  If high-VAF calls are hairpin-depleted, the low band MUST look depleted. That is a")
print("  property of the calling rule, not of the estimator.")
