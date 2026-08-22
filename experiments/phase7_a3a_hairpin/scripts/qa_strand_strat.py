#!/usr/bin/env python
"""STRAND SYMMETRY of the Y130G signal -- motivated directly by BUG 4.

Bug 4 was `s6_pileup` counting only uppercase T, which drops half the reads
STRAND-ASYMMETRICALLY. That bug is fixed, but its existence is the reason to ask whether the
headline signal is strand-balanced: a residual strand-specific counting or alignment defect
would show up as an asymmetry between the two strands, while a real deaminase preference for
hairpin loops should appear on both.

The universe carries a per-site `strand` field (0 = C on plus, 1 = C on minus). Split the
same estimator by it. Both strands should show the effect; a large asymmetry is a flag.

Also reports the ALT-READ strand balance within the editor-specific sites themselves
(alt_fwd vs alt_rev), which is the direct descendant of bug 4's failure mode.
"""
import numpy as np, sys, os
FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED, CAL, MASK = sys.argv[1], sys.argv[2], "Parent"
STEM = 6

acc = {0: dict(nel=0, nhp=0, e=0, ehp=0, c=0, chp=0, af=0, ar=0),
       1: dict(nel=0, nhp=0, e=0, ehp=0, c=0, chp=0, af=0, ar=0)}
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    st, hp = U["strand"], U["stem"] >= STEM
    elig = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["cov"] >= 15) & (P["alt"] == 0)
    espec = elig & (E["alt"] >= 2) & (K["alt"] == 0)
    cspec = elig & (K["alt"] >= 2) & (E["alt"] == 0)
    for s in (0, 1):
        m = st == s
        a = acc[s]
        a["nel"] += int((elig & m).sum()); a["nhp"] += int((elig & m & hp).sum())
        a["e"]   += int((espec & m).sum()); a["ehp"] += int((espec & m & hp).sum())
        a["c"]   += int((cspec & m).sum()); a["chp"] += int((cspec & m & hp).sum())
        a["af"]  += int(E["alt_fwd"][espec & m].sum()); a["ar"] += int(E["alt_rev"][espec & m].sum())
    del U, E, K, P

print(f"=== {ED} vs {CAL} -- STRAND-STRATIFIED (stem>={STEM}) ===")
print(f"  {'strand':>8} {'n_elig':>13} {'p_bg':>8} {'editor':>8} {'n_hp':>7} {'calib':>8} {'n_hp':>7} {'ed-cal':>8}")
vals = {}
for s in (0, 1):
    a = acc[s]
    pbg = a["nhp"] / max(a["nel"], 1)
    fe = (a["ehp"] / max(a["e"], 1)) / pbg
    fc = (a["chp"] / max(a["c"], 1)) / pbg
    vals[s] = (fe, fc)
    lab = "0 (plus)" if s == 0 else "1 (minus)"
    print(f"  {lab:>8} {a['nel']:13,} {pbg:8.5f} {fe:8.3f} {a['ehp']:7,} {fc:8.3f} {a['chp']:7,} {fe-fc:+8.3f}")
d0 = vals[0][0] - vals[0][1]; d1 = vals[1][0] - vals[1][1]
print(f"\n  ed-cal on plus {d0:+.3f}   on minus {d1:+.3f}   asymmetry {abs(d0-d1):.3f}")
print("  Both strands should carry it. A large asymmetry is the signature bug 4 left behind.")
print("\n=== ALT-READ balance inside the editor-specific sites (bug 4's own failure mode) ===")
for s in (0, 1):
    a = acc[s]
    print(f"  strand {s}: alt_fwd {a['af']:,}  alt_rev {a['ar']:,}  ratio {a['af']/max(a['ar'],1):.4f}")
