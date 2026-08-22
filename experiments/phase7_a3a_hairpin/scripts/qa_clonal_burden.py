#!/usr/bin/env python
"""QA of my OWN clonal-burden numbers from last tick.

I reported ratios like "A3A-Y130F-clone2 vs nCas9-clone1 = 1.08x, z=1.38" from RAW COUNTS of
clonal variants. That violates the standing rule: the burden bar applies coverage-MATCHED.
Raw counts compare two samples over different eligible sets at different depths, which is the
exact confound Endpoint A's banding exists to remove -- I applied correct statistics to the
wrong quantity.

Recomputing burden as a RATE over each sample's OWN eligible set, within matched coverage
bands, restricted to CLONAL variants (VAF >= 0.35), with the cov-ratio audit printed so a band
that did not equalise depth can be discarded rather than read.
"""
import sys
import numpy as np

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
BANDS = [(8, 15), (15, 25), (25, 35), (35, 60)]
VAF_MIN = 0.35
ED = sys.argv[1] if len(sys.argv) > 1 else "P66-A3A-Y130F-clone2"
CAL = sys.argv[2] if len(sys.argv) > 2 else "nCas9-clone1"

acc = {b: [0, 0, 0, 0, 0, 0] for b in BANDS}   # ed_hit, ed_elig, cal_hit, cal_elig, ed_cov, cal_cov
for c in CH:
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15)
    for nm, D, o in (("ed", E, 0), ("cal", K, 2)):
        cov, alt = D["cov"], D["alt"]
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(cov > 0, alt / np.maximum(cov, 1), 0.0)
        for lo, hi in BANDS:
            b = sil & (cov >= lo) & (cov < hi)
            hit = b & (alt >= 2) & (v >= VAF_MIN)
            a = acc[(lo, hi)]
            a[o] += int(hit.sum()); a[o + 1] += int(b.sum())
            a[4 if o == 0 else 5] += int(cov[b].sum())

print(f"=== CLONAL burden (VAF>={VAF_MIN}), rate per Mb of each sample's OWN eligible set ===")
print(f"    editor {ED}   calibrator {CAL}")
print(f"  {'band':>14} {'ed/Mb':>8} {'cal/Mb':>8} {'ratio':>7} {'ed cov':>7} {'cal cov':>8} "
      f"{'cov ratio':>10} {'n_ed':>7}")
for lo, hi in BANDS:
    eh, ee, kh, ke, ec, kc = acc[(lo, hi)]
    if ee < 100000 or ke < 100000:
        print(f"  {f'({lo},{hi})':>14}  -- too few eligible sites --"); continue
    er = eh / ee * 1e6; kr = kh / ke * 1e6
    cr = (ec / ee) / (kc / ke)
    flag = "" if 0.98 <= cr <= 1.02 else "   <-- UNMATCHED DEPTH, bar does not apply"
    print(f"  {f'({lo},{hi})':>14} {er:>8.1f} {kr:>8.1f} {er/kr if kr else float('nan'):>6.2f}x "
          f"{ec/ee:>7.2f} {kc/ke:>8.2f} {cr:>10.3f} {eh:>7,}{flag}")
print(f"\n  Compare with the RAW-COUNT ratio I quoted last tick for this pair.")
print(f"  A band's ratio is only readable when its cov ratio sits in 0.98-1.02.")
