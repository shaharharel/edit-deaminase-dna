#!/usr/bin/env python
"""GATE A0 v3: the noise-floor criteria were DEPTH-CONFOUNDED. Fix the generator.

MEASURED across the seven samples with 23/23 counts:
    corr(coverage, %<VAF.05) = 0.9255
At every depth in this set, alt=1 sits below VAF 0.05 and alt=2 above it, so
"%<VAF.05" is essentially "fraction of alt>=1 sites that are singletons" -- and
singletons scale with depth x error rate while true variants do not. That criterion
is a DEPTH DETECTOR, not a quality detector, and it must not stand alone.

The COUNT criterion survives, but only once depth-adjusted. Fitting log10(alt>=1) on
coverage using ONLY the five samples the gate calls healthy (r = 0.9493):
    nCas9-clone2      0.97x     A3A-Y130F-clone5  0.95x
    nCas9-clone1      1.07x     Parent            1.02x
    A3A-Y130F-clone2  0.99x
    ---------------------------------------------
    D10A-clone1       5.15x  <-- disqualified
    D10A-clone6       2.96x  <-- disqualified
Sound samples land in 0.95-1.07x. The D10A clones sit multiples above. THAT is the
real signal, and the raw peer-median ratio only found it by luck, because the D10A
clones happen also to be the deepest.

v3 therefore: (a) reports the depth-adjusted residual as the PRIMARY criterion,
(b) keeps %<VAF.05 but labels it depth-confounded and advisory, (c) still refuses.
"""
import os, sys, shutil

BASE = os.environ.get("A3A_BASE", "/mnt/data/a3a")
p = f"{BASE}/qualify_calibrator.py"
src = open(p).read()
SENT = "# --- v3 depth-adjusted residual (added 2026-08-22) ---"
if SENT in src:
    print("  already v3"); sys.exit(0)

ANCHOR = "if f_cal > MAX_SUBVAF:"
if ANCHOR not in src:
    print("  *** anchor not found -- NOT patched ***"); sys.exit(1)

NEW = SENT + '''
# The %<VAF.05 criterion below correlates 0.9255 with coverage across the seven samples
# measured, because at these depths alt=1 is under VAF 0.05 and alt=2 is over it, so the
# statistic is "fraction of alt>=1 that are singletons" -- which scales with depth. It is
# kept as ADVISORY. The primary criterion is now the depth-adjusted residual: fit
# log10(alt>=1) on coverage over the PEERS ONLY, and ask how far the calibrator sits above
# its own expectation. Sound samples land at 0.95-1.07x; the two D10A clones at 2.96x and
# 5.15x.
MAX_RESID = 2.0
resid = None
if len(peers) >= 3:
    import numpy as _np
    _c = _np.array([v[2] for v in peers.values()], dtype=float)
    _n = _np.array([v[0] for v in peers.values()], dtype=float)
    _keep = (_n > 0) & (_c > 0)
    if _keep.sum() >= 3:
        _b, _a = _np.polyfit(_c[_keep], _np.log10(_n[_keep]), 1)
        _exp = 10 ** (_a + _b * d_cal)
        resid = n_cal / _exp
        print(f"  depth-adjusted: expected {_exp:,.0f} alt>=1 at cov {d_cal:.2f}, "
              f"observed {n_cal:,} -> {resid:.2f}x above the peer depth trend "
              f"(limit {MAX_RESID:.1f}x)")
    else:
        print("  depth-adjusted: too few usable peers to fit a trend -- falling back")
else:
    print("  depth-adjusted: fewer than 3 peers -- cannot fit a depth trend")

if resid is not None and resid > MAX_RESID:
    fail.append(f"{resid:.2f}x above the peer DEPTH TREND (limit {MAX_RESID:.1f}x) "
                f"-- this is the depth-adjusted criterion, not the raw ratio")

''' + ANCHOR
src = src.replace(ANCHOR, NEW, 1)
src = src.replace(
    'fail.append(f"{100*f_cal:.1f}% of alt>=1 sites below VAF 0.05 (limit {100*MAX_SUBVAF:.0f}%)")',
    'fail.append(f"{100*f_cal:.1f}% of alt>=1 sites below VAF 0.05 (limit {100*MAX_SUBVAF:.0f}%) '
    '[ADVISORY: this statistic correlates 0.93 with coverage]")', 1)
shutil.copy(p, p + ".prev3")
open(p, "w").write(src)
compile(open(p).read(), p, "exec")
print(f"  GATE FIXED to v3 (backup .prev3), compiles")
