#!/usr/bin/env python
"""CHECK 4 applied to the EDITOR arms, which it never has been.

GC-decile and complexity stratification have been run on the PCAWG model all night. The editor
endpoint -- hairpin enrichment among Parent-masked specific sites -- has never been stratified
at all, and it is just as exposed: hairpins need inverted repeats, GC-rich sequence forms more
stable stems, and variant calling itself has GC-dependent sensitivity. If the enrichment lives
in one GC stratum it is a composition effect, not a hairpin effect.

Runs the editor and its deaminase-free calibrator side by side in every decile, so a shared
GC artefact shows up as BOTH arms moving together while a real difference does not.
"""
import sys
import numpy as np

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED = sys.argv[1] if len(sys.argv) > 1 else "P66-A3A-Y130F-clone2"
CAL = sys.argv[2] if len(sys.argv) > 2 else "nCas9-clone1"
STEM = 6          # stem>=8 has too few sites to stratify; check 6 forbids pretending otherwise

gc_all, hp_all, sp_ed, sp_cal = [], [], [], []
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    sil = (P["alt"] == 0) & (P["cov"] >= 15) & (E["cov"] >= 8) & (K["cov"] >= 8)
    idx = np.flatnonzero(sil)
    gc_all.append(U["local_gc"][idx].astype(np.float32))
    hp_all.append((U["stem"][idx] >= STEM))
    sp_ed.append(E["alt"][idx] >= 2)
    sp_cal.append(K["alt"][idx] >= 2)
gc = np.concatenate(gc_all); hp = np.concatenate(hp_all)
se = np.concatenate(sp_ed); sk = np.concatenate(sp_cal)
print(f"jointly eligible sites: {len(gc):,}   stem>={STEM} background {hp.mean():.5f}")
print(f"editor {ED} specific {int(se.sum()):,}   calibrator {CAL} specific {int(sk.sum()):,}")

qs = np.unique(np.percentile(gc, np.arange(0, 101, 10)))
print(f"\n  {'decile':>6} {'GC range':>13} {'n_elig':>11} {'p_bg':>8} "
      f"{'editor':>9} {'n_hp':>7} {'calib':>9} {'n_hp':>7}")
for i in range(len(qs) - 1):
    lo, hi = qs[i], qs[i + 1]
    m = (gc >= lo) & (gc < hi) if i < len(qs) - 2 else (gc >= lo) & (gc <= hi)
    if m.sum() < 100000:
        continue
    pbg = hp[m].mean()
    ee = (hp[m & se].mean() / pbg) if (m & se).sum() > 50 and pbg > 0 else float("nan")
    ek = (hp[m & sk].mean() / pbg) if (m & sk).sum() > 50 and pbg > 0 else float("nan")
    print(f"  {i:>6} {lo:5.3f}-{hi:<7.3f} {int(m.sum()):>11,} {pbg:>8.5f} "
          f"{ee:>8.3f}x {int((m & se & hp).sum()):>7,} {ek:>8.3f}x {int((m & sk & hp).sum()):>7,}")
print("\n  A shared GC artefact moves BOTH columns together. A real editor difference does not.")
print("  n_hp is printed per decile because check 6 applies here as much as anywhere.")

# p_bg spans 2.9x across GC deciles, so ANY shift in the GC composition of the called sites
# moves the pooled enrichment directly. Same shape as the depth confound that banding fixed;
# same remedy -- stratify, then pool by Mantel-Haenszel.
print("\n=== GC-ADJUSTED pooling (Mantel-Haenszel over the deciles above) ===")
for lab, sel in (("editor", se), ("calibrator", sk)):
    num = den = 0.0
    for i in range(len(qs) - 1):
        lo, hi = qs[i], qs[i + 1]
        m = (gc >= lo) & (gc < hi) if i < len(qs) - 2 else (gc >= lo) & (gc <= hi)
        if m.sum() < 100000:
            continue
        A = float((m & sel & hp).sum()); B = float((m & sel & ~hp).sum())
        C = float((m & ~sel & hp).sum()); D = float((m & ~sel & ~hp).sum())
        T = A + B + C + D
        if T <= 0: continue
        num += A * D / T; den += B * C / T
    crude_h = float((sel & hp).sum()); crude_f = float((sel & ~hp).sum())
    bg_h = float((~sel & hp).sum()); bg_f = float((~sel & ~hp).sum())
    crude = (crude_h / crude_f) / (bg_h / bg_f)
    print(f"  {lab:11s} crude OR {crude:.3f}   GC-adjusted MH {num/den:.3f}   "
          f"shift {100*(num/den - crude)/crude:+.1f}%")
print("  If adjustment moves the editor materially and the calibrator little, the pooled")
print("  editor number was carrying a GC-composition effect.")
