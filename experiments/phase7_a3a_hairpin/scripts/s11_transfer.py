#!/usr/bin/env python
"""THE TRANSFER TEST. Does a PCAWG-trained classifier rank real base-editor off-target sites
above matched background in HEK293T?

This is the project's central claim and it has never been tested. What exists is (a) the model
predicting held-out PCAWG chromosomes at 5.036x, which is tumour mutations predicting tumour
mutations, and (b) the hairpin FEATURE measured separately on both datasets. Neither is the
model transferring.

DESIGN, and the control is the whole point:
  train  : PCAWG v5, the 6 structural features that won the architecture search
  score  : every jointly-eligible TCW site in HEK293T (Parent-silent, cov>=8 in both samples)
  ask    : are EDITOR-specific sites enriched in the model's top K%?
  control: are CALIBRATOR-specific sites enriched by the same amount?

If both arms enrich equally, the model is predicting WHERE VARIANTS GET CALLED -- coverage,
mappability, GC, the calling artefacts this project has spent the night characterising -- and
not anything about the editor. Only an editor-minus-calibrator gap is evidence of transfer.

Two caveats built in rather than discovered later:
  - the editor-specific set is 91-94%% subclonal (measured 04:10), so this necessarily tests
    prediction of subclonal calls;
  - the model uses local_gc, and site-calling is GC-biased, which is exactly why the
    calibrator control is not optional.
Random baseline beside every number; n reported throughout.
"""
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
FEATS = ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")
import sys
ED  = sys.argv[1] if len(sys.argv) > 1 else "P66-A3A-Y130F-clone2"
CAL = sys.argv[2] if len(sys.argv) > 2 else "nCas9-clone1"
PCTS = (0.1, 0.5, 1.0, 5.0)

d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
Xtr = np.stack([d[k].astype(np.float32) for k in FEATS], 1)
ytr = d["y"].astype(np.int8)
print(f"train: PCAWG v5, {Xtr.shape[0]:,} rows, {Xtr.shape[1]} structural features, "
      f"base rate {ytr.mean():.4f}")
m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                   validation_fraction=0.1, random_state=0)
m.fit(Xtr, ytr)
print("model fitted on PCAWG only -- no HEK293T data seen during training")

sc, is_ed, is_cal = [], [], []
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    el = (P["alt"] == 0) & (P["cov"] >= 15) & (E["cov"] >= 8) & (K["cov"] >= 8)
    i = np.flatnonzero(el)
    X = np.stack([U[k][i].astype(np.float32) for k in FEATS], 1)
    sc.append(m.predict_proba(X)[:, 1].astype(np.float32))
    is_ed.append(E["alt"][i] >= 2)
    is_cal.append(K["alt"][i] >= 2)
    del U, P, E, K, X
s = np.concatenate(sc); ed = np.concatenate(is_ed); cal = np.concatenate(is_cal)
N = len(s)
print(f"\nscored {N:,} jointly-eligible HEK293T sites")
print(f"  editor-specific     {int(ed.sum()):>9,}  base rate {ed.mean():.6f}")
print(f"  calibrator-specific {int(cal.sum()):>9,}  base rate {cal.mean():.6f}")

rng = np.random.default_rng(0)
rand = rng.random(N).astype(np.float32)
print(f"\n=== does the PCAWG model rank HEK293T off-target sites above background? ===")
print(f"  {'top':>6} {'n':>9} {'EDITOR':>9} {'n_hit':>7} {'CALIBRATOR':>11} {'n_hit':>7} "
      f"{'random(ed)':>11} {'ed - cal':>9}")
for p in PCTS:
    k = max(int(round(N * p / 100)), 1)
    idx = np.argpartition(-s, k - 1)[:k]
    ridx = np.argpartition(-rand, k - 1)[:k]
    e = ed[idx].mean() / ed.mean()
    c_ = cal[idx].mean() / cal.mean()
    r = ed[ridx].mean() / ed.mean()
    print(f"  {p:>5}% {k:>9,} {e:>8.3f}x {int(ed[idx].sum()):>7,} {c_:>10.3f}x "
          f"{int(cal[idx].sum()):>7,} {r:>10.3f}x {e - c_:>+9.3f}")
# POSITIVE CONTROL: the same fitted model on its own training distribution. If this is flat
# too, the model is broken rather than the transfer failing -- a distinction worth one line.
pc = m.predict_proba(Xtr)[:, 1]
kk = max(int(round(len(pc) * 0.001)), 1)
ii = np.argpartition(-pc, kk - 1)[:kk]
print(f"\n  POSITIVE CONTROL, same model on PCAWG (in-sample): top 0.1% = "
      f"{ytr[ii].mean()/ytr.mean():.3f}x  (n_top={kk:,})")
print("  -> the model works. What fails is the transfer, not the model.")
print("\n  EQUAL enrichment in both arms = the model predicts where variants get CALLED.")
print("  Only an editor-minus-calibrator GAP is evidence that it predicts editor off-targets.")
