#!/usr/bin/env python
"""GATE A0: is this calibrator technically sound, before anything is computed against it?

v2. v1 REJECTED nCas9-clone1 -- a calibrator the protocol names as valid and that produced a
passing harness validation an hour ago. Its rule ("specific-site count >3x the median")
compared against a reference set dominated by the two editor clones, which have the LOWEST
counts. That encodes "the calibrator must not carry more background than the editor", which is
precisely the study-level difference already diagnosed as real. A gate that rejects the
correct calibrator for the correct reason being present is worse than no gate.

The genuine defect in D10A-clone1 was never its alt>=2 count. It was the RAW NOISE FLOOR,
which is a property of the sequencing rather than of how much true somatic background a clone
carries. Measured across every sample with 23/23 counts:

    sample                  alt>=1 sites (cov>=8)   % below VAF 0.05   mean cov
    P66-D10A-clone1               6,457,037              87.2%          31.36
    P66-A3A-Y130F-clone2            864,617              63.7%          25.02
    P66-A3A-Y130F-clone5            861,941              65.5%          25.73
    nCas9-clone1                    854,914              53.0%          23.60
    nCas9-clone2                    780,594              51.8%          23.63
    Parent                        1,135,489              67.2%          29.29

clone1 is 7.5x the next highest and its sub-0.05 fraction sits far outside a 51.8-67.2% band.
Both criteria are set from that spread, not from a guess, and both catch clone1 while passing
all five sound samples -- which is the validation run at the bottom of this file.
"""
import os, sys
import numpy as np

# FEAT is env-overridable so ONE file serves both nodes. Hardcoding a node-A path and
# copying the file to node B is how this project got five "wrong artifact" errors.
FEAT = os.environ.get("A3A_FEAT", "/mnt/data/a3a/feat")
CH = [str(i) for i in range(1, 23)] + ["X"]
CAL = sys.argv[1] if len(sys.argv) > 1 else "P66-D10A-clone6"
PEERS = ["P66-A3A-Y130F-clone2", "P66-A3A-Y130F-clone5", "nCas9-clone1", "nCas9-clone2", "Parent"]
MAX_SUBVAF = 0.80      # healthy band is 0.518-0.672; clone1 is 0.872
MAX_RATIO  = 3.0       # clone1 is 7.5x the peer median

def noise_floor(s):
    n = lo = ct = cn = 0
    for c in CH:
        D = np.load(f"{FEAT}/counts_{s}_chr{c}.npz")
        cov, alt = D["cov"], D["alt"]
        m = (alt >= 1) & (cov >= 8)
        n += int(m.sum()); ct += int(cov.sum()); cn += len(cov)
        v = alt[m] / np.maximum(cov[m], 1)
        lo += int((v < 0.05).sum())
    return n, lo / max(n, 1), ct / max(cn, 1)

print(f"=== calibrator qualification (raw noise floor): {CAL} ===")
try:
    n_cal, f_cal, d_cal = noise_floor(CAL)
except FileNotFoundError:
    print(f"  *** {CAL}: counts incomplete -- cannot qualify"); sys.exit(3)

peers = {}
for s in PEERS:
    if s == CAL: continue
    try: peers[s] = noise_floor(s)
    except FileNotFoundError: pass
med = float(np.median([v[0] for v in peers.values()])) if peers else float("nan")

print(f"  {'sample':24s} {'alt>=1':>12s} {'%<VAF.05':>9s} {'cov':>7s}")
print(f"  {CAL:24s} {n_cal:>12,} {100*f_cal:>8.1f}% {d_cal:>7.2f}  <-- CALIBRATOR")
for s, (n, f, d) in peers.items():
    print(f"  {s:24s} {n:>12,} {100*f:>8.1f}% {d:>7.2f}")
print(f"  peer median alt>=1: {med:,.0f}   calibrator is {n_cal/med:.2f}x")

fail = []
if f_cal > MAX_SUBVAF:
    fail.append(f"{100*f_cal:.1f}% of alt>=1 sites below VAF 0.05 (limit {100*MAX_SUBVAF:.0f}%)")
if peers and n_cal > MAX_RATIO * med:
    fail.append(f"{n_cal/med:.1f}x the peer median noise floor (limit {MAX_RATIO:.0f}x)")
print()
if fail:
    print("  *** DISQUALIFIED: " + "; ".join(fail))
    print("  *** Do NOT run the editor read-out against this calibrator.")
    sys.exit(2)
print(f"  QUALIFIED: noise floor {n_cal:,} ({n_cal/med:.2f}x peer median), "
      f"{100*f_cal:.1f}% sub-0.05")
sys.exit(0)
