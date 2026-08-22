#!/usr/bin/env python
"""THE MISSING CONTROL: the editor showed a mirror-image alt-read skew tied to the C's strand
(1.0903 on plus-strand C, 0.8880 on minus-strand C). That is bug 4's exact family.

But I measured it ONLY on the editor. If the CALIBRATOR shows the same mirror skew, it is a
property of the pipeline or of the reference base, harmless and shared. If ONLY the editor
shows it, it is a real concern for the editor-specific set.

Measure alt_fwd/alt_rev, split by the C's genomic strand, for:
    editor-specific sites, calibrator-specific sites, and ALL eligible called sites in each
So the comparison is like-for-like.
"""
import numpy as np, sys
FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED, CAL, MASK = sys.argv[1], sys.argv[2], "Parent"

acc = {}
for key in ("ed_spec", "cal_spec", "ed_all", "cal_all"):
    acc[key] = {0: [0, 0], 1: [0, 0]}
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    st = U["strand"]
    elig = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["cov"] >= 15) & (P["alt"] == 0)
    sets = {"ed_spec":  (elig & (E["alt"] >= 2) & (K["alt"] == 0), E),
            "cal_spec": (elig & (K["alt"] >= 2) & (E["alt"] == 0), K),
            "ed_all":   (elig & (E["alt"] >= 2), E),
            "cal_all":  (elig & (K["alt"] >= 2), K)}
    for key, (m, D) in sets.items():
        for s in (0, 1):
            mm = m & (st == s)
            acc[key][s][0] += int(D["alt_fwd"][mm].sum())
            acc[key][s][1] += int(D["alt_rev"][mm].sum())
    del U, E, K, P

print(f"=== alt_fwd / alt_rev, split by the C's genomic strand ===")
print(f"  {'set':12s} {'C on plus':>26} {'C on minus':>26}")
print(f"  {'':12s} {'fwd':>10} {'rev':>9} {'ratio':>6} {'fwd':>10} {'rev':>9} {'ratio':>6}")
for key in ("ed_spec", "cal_spec", "ed_all", "cal_all"):
    a = acc[key]
    r0 = a[0][0] / max(a[0][1], 1); r1 = a[1][0] / max(a[1][1], 1)
    print(f"  {key:12s} {a[0][0]:10,} {a[0][1]:9,} {r0:6.3f} {a[1][0]:10,} {a[1][1]:9,} {r1:6.3f}")
e0 = acc["ed_spec"][0][0]/max(acc["ed_spec"][0][1],1); e1 = acc["ed_spec"][1][0]/max(acc["ed_spec"][1][1],1)
c0 = acc["cal_spec"][0][0]/max(acc["cal_spec"][0][1],1); c1 = acc["cal_spec"][1][0]/max(acc["cal_spec"][1][1],1)
print(f"\n  editor      mirror skew: {e0:.3f} / {e1:.3f}   spread {abs(e0-e1):.3f}")
print(f"  calibrator  mirror skew: {c0:.3f} / {c1:.3f}   spread {abs(c0-c1):.3f}")
print("\n  If the calibrator shows the SAME mirror pattern, it is a shared pipeline property")
print("  and harmless. If only the editor does, it is a real concern for the editor set.")
