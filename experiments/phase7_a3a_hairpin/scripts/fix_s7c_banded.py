#!/usr/bin/env python
"""s7c: ENDPOINT B WITHIN MATCHED COVERAGE BANDS -- the fix that replaces the VAF floor.

WHY THE VAF FLOOR WAS NOT THE ANSWER. alt>=2 is a depth-dependent VAF threshold (VAF>=0.05
at 40x, >=0.087 at 23x), so arms of unequal depth call different site populations. I patched
in A3A_VAF_MIN, then ran the nCas9-vs-nCas9 validation where both clones are deaminase-free
and must be indistinguishable:
    stem>=6   no floor  1.102 / 1.051   (n_hp 1,480)
              VAF>=0.15 1.438 / 1.159   (n_hp   165)
The floor made two identical-treatment clones diverge and cost 90% of the power. It traded a
depth confound for a clone-VAF confound.

WHAT THIS DOES INSTEAD. Endpoint A already solves this problem for burden: compare within
coverage bands, and audit that the band actually equalised depth. Endpoint B gets the same
treatment. Within each band the eligible set requires BOTH samples to have coverage in that
band, so calling sensitivity is equalised at the same sites for both arms before any hairpin
fraction is taken. Per-band enrichments are reported with n, and pooled by Mantel-Haenszel
-- the stratified summary a band table implies but does not state.

This is the same instrument that cleared the PCAWG coverage confound (MH 1.323 vs crude
1.349 at stem>=8), applied to the editor arms.

Output is ADDITIVE: the existing pooled endpoint B stays, so every prior number remains
reproducible and the banded version sits beside it for comparison.
"""
import sys, os, shutil

BASE = sys.argv[1] if len(sys.argv) > 1 else "/mnt/data/a3a"
p = f"{BASE}/s7c_editor.py"
src = open(p).read()

SENT = "# --- banded endpoint B (added 2026-08-21) ---"
if SENT in src:
    print("  already patched"); sys.exit(0)

# 1. accumulator, declared beside the existing ones
A1 = "burden = {b: [0, 0, 0, 0] for b in BANDS}"
if A1 not in src:
    print("  *** burden accumulator not found - NOT PATCHED ***"); sys.exit(1)
src = src.replace(A1, A1 + "\n" + SENT + """
# [band][stem][w] = [n_hp_spec, n_spec, n_hp_elig, n_elig]
hp_band = {b: {s: {"ed": [0, 0, 0, 0], "cal": [0, 0, 0, 0]} for s in STEMS} for b in BANDS}""", 1)

# 2. fill it in the per-chromosome loop, right after the pooled endpoint B block
A2 = """    for w, D in (("ed", E), ("cal", K)):
        elig = silent & (D["cov"] >= 8)
        spec = elig & _called(D, 0)"""
if A2 not in src:
    print("  *** endpoint B loop not found - NOT PATCHED ***"); sys.exit(1)
NEW2 = """    # banded endpoint B: eligibility requires BOTH samples inside the same coverage
    # band, so site-calling sensitivity is equalised before any hairpin fraction is taken.
    for lo, hi in BANDS:
        both = silent & (E["cov"] >= lo) & (E["cov"] < hi) & (K["cov"] >= lo) & (K["cov"] < hi)
        if not both.any():
            continue
        for w, D in (("ed", E), ("cal", K)):
            sp = both & _called(D, 0)
            for s in STEMS:
                h = stem >= s
                a = hp_band[(lo, hi)][s][w]
                a[0] += int((sp & h).sum()); a[1] += int(sp.sum())
                a[2] += int((both & h).sum()); a[3] += int(both.sum())

""" + A2
src = src.replace(A2, NEW2, 1)

# 3. report it after the existing endpoint B table
A3 = 'print("\\nENDPOINT B by trinucleotide context'
if A3 not in src:
    print("  *** context table anchor not found - NOT PATCHED ***"); sys.exit(1)
NEW3 = '''print("ENDPOINT B -- BANDED (depth-matched; alt>=2 is a depth-dependent VAF cut)")
print("   stem   band          editor    calib   n_hp_ed  n_hp_cal   p_bg")
for s in STEMS:
    _num_e = _den_e = _num_c = _den_c = 0.0
    for lo, hi in BANDS:
        a = hp_band[(lo, hi)][s]
        ed, cal = a["ed"], a["cal"]
        if ed[3] < 500 or ed[1] < 20 or cal[1] < 20:
            continue
        pbg = ed[2] / ed[3]
        e_ed = (ed[0] / ed[1]) / pbg if ed[1] else float("nan")
        e_cal = (cal[0] / cal[1]) / pbg if cal[1] else float("nan")
        print(f"   {s:>4}   ({lo},{hi})".ljust(20)
              + f"{e_ed:8.3f} {e_cal:8.3f} {ed[0]:>9,} {cal[0]:>9,}  {pbg:.5f}")
        # Mantel-Haenszel across bands, hairpin vs flat within each arm
        for _n, _d, arr in ((0, 0, ed), (1, 1, cal)):
            A = arr[0]; B = arr[1] - arr[0]
            C = arr[2] - arr[0]; Dd = arr[3] - arr[1] - C
            T = A + B + C + Dd
            if T <= 0: continue
            if arr is ed: _num_e += A * Dd / T; _den_e += B * C / T
            else: _num_c += A * Dd / T; _den_c += B * C / T
    if _den_e > 0 and _den_c > 0:
        print(f"   {s:>4}   MANTEL-HAENSZEL pooled OR   editor {_num_e/_den_e:.3f}"
              f"   calibrator {_num_c/_den_c:.3f}")
print("   (bands with <500 eligible or <20 called sites are omitted, not silently pooled)")
print()
''' + A3
src = src.replace(A3, NEW3, 1)

shutil.copy(p, p + ".preband")
open(p, "w").write(src)
compile(open(p).read(), p, "exec")
print("  s7c_editor.py PATCHED: banded endpoint B added (backup .preband)")
print("  pooled endpoint B is unchanged, so every prior number stays reproducible")
