#!/usr/bin/env python
"""QA of my OWN transfer result. A strong negative deserves the scrutiny a positive would get.

THE ALTERNATIVE EXPLANATION I have not excluded: the model ranks hairpin-rich sites high,
hairpins are commoner in AT-rich sequence (measured tonight: p_bg spans 2.9x across GC
deciles), and AT-rich sequence has lower mappability and coverage. A variant needs COVERAGE to
be called at all. So the model's top-scoring sites may be systematically the ones where a
variant is HARDEST to call -- which would produce apparent depletion of called sites with no
involvement of the editor whatsoever.

If that is what happened, my "the transfer fails" conclusion is confounded in exactly the way
everything else tonight was, and the test needs coverage-matching before it means anything.

Checks: coverage and GC of the model's top-ranked sites vs the whole eligible set, and the
transfer result recomputed WITHIN coverage bands.
"""
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
FEATS = ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")
ED, CAL = "P66-A3A-Y130F-clone2", "nCas9-clone1"
BANDS = [(8, 15), (15, 25), (25, 35), (35, 60)]

d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                   validation_fraction=0.1, random_state=0)
m.fit(np.stack([d[k].astype(np.float32) for k in FEATS], 1), d["y"].astype(np.int8))

sc, ed_l, cal_l, cov_l, gc_l, tri_l = [], [], [], [], [], []
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    el = (P["alt"] == 0) & (P["cov"] >= 15) & (E["cov"] >= 8) & (K["cov"] >= 8)
    i = np.flatnonzero(el)
    sc.append(m.predict_proba(np.stack([U[k][i].astype(np.float32) for k in FEATS], 1))[:, 1].astype(np.float32))
    ed_l.append(E["alt"][i] >= 2); cal_l.append(K["alt"][i] >= 2)
    cov_l.append(E["cov"][i].astype(np.float32)); gc_l.append(U["local_gc"][i].astype(np.float32))
    tri_l.append(U["tri"][i])
    del U, P, E, K
s = np.concatenate(sc); ed = np.concatenate(ed_l); cal = np.concatenate(cal_l)
cov = np.concatenate(cov_l); gc = np.concatenate(gc_l); tri = np.concatenate(tri_l)
N = len(s)

print(f"=== 0. is the scoring set the same UNIVERSE the model was trained on? ===")
u, ct = np.unique(tri[:2000000], return_counts=True)
print(f"  trinucleotides in the scoring set: {dict(zip(u.tolist(), (ct/ct.sum()).round(4).tolist()))}")
print(f"  (the model was trained on TCA/TCT only -- anything else here is out of distribution)")

print(f"\n=== 1. does the model rank LOW-COVERAGE sites high? (the confound) ===")
print(f"  {'set':>18} {'n':>12} {'mean cov':>9} {'mean GC':>9}")
print(f"  {'all eligible':>18} {N:>12,} {cov.mean():>9.2f} {gc.mean():>9.4f}")
for p in (0.1, 1.0, 5.0):
    k = max(int(round(N * p / 100)), 1)
    idx = np.argpartition(-s, k - 1)[:k]
    print(f"  {f'model top {p}%':>18} {k:>12,} {cov[idx].mean():>9.2f} {gc[idx].mean():>9.4f}")

print(f"\n=== 2. the transfer test WITHIN coverage bands (removes the confound if it exists) ===")
print(f"  {'band':>12} {'n_elig':>12} {'editor':>9} {'n_hit':>7} {'calib':>9} {'n_hit':>7}")
for lo, hi in BANDS:
    b = (cov >= lo) & (cov < hi)
    if b.sum() < 500000: continue
    sb = s[b]; edb = ed[b]; calb = cal[b]
    k = max(int(round(b.sum() * 0.01)), 1)          # top 1% within the band
    idx = np.argpartition(-sb, k - 1)[:k]
    e = edb[idx].mean() / edb.mean() if edb.mean() > 0 else float("nan")
    c_ = calb[idx].mean() / calb.mean() if calb.mean() > 0 else float("nan")
    print(f"  {f'({lo},{hi})':>12} {int(b.sum()):>12,} {e:>8.3f}x {int(edb[idx].sum()):>7,} "
          f"{c_:>8.3f}x {int(calb[idx].sum()):>7,}")
print("  (top 1% WITHIN each band, so coverage cannot drive the comparison)")
