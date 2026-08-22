#!/usr/bin/env python
"""TOP-K PANEL CAPTURE, PCAWG-trained model -> Y130G. The test that does not exist yet.

WHY THIS AND NOT A REPEAT OF THE TRANSFER TEST. The panel metric IS the enrichment metric:
a panel holding K% of sites that captures K% of the editor's mutations scores 1.00x and has
no skill. The transfer test already measured it -- 0.954x and 1.001x at top 5% -- so for
A3A-Y130F the answer is in and it is no.

But that test was only ever run against A3A-Y130F, THE ARM WITH NO HAIRPIN SIGNAL. A model
whose learned signal is hairpin geometry cannot rank sites in a sample that has no hairpin
excess. It can be expected to rank them in one that does -- and Y130G shows +0.318/+0.396
GC-adjusted MH, replicated over two clones and surviving GC, complexity, strand, mix-shift
and read-orientation audits.

So: train on PCAWG only, score every jointly-eligible Y130G site, and ask what fraction of
the editor's mutation burden a top-K panel captures.

THE CALIBRATOR COLUMN IS THE WHOLE TEST. nCas9-clone1 is deaminase-free and burden-matched to
Y130G within 1%. If the panel captures the editor and the calibrator EQUALLY, the model is
finding where variants get CALLED, not where the editor acts. Only an editor-minus-calibrator
gap is transfer.

Pre-registered read, written before the numbers exist:
  capture ratio ~1.0 in both columns          -> no transfer, same answer as A3A-Y130F
  editor > calibrator at small K              -> the PCAWG model carries usable hairpin
                                                 knowledge into an arm that has hairpin signal
  editor < calibrator                         -> anti-transfer, as seen for A3A-Y130F
"""
import numpy as np, os, sys, time
from sklearn.ensemble import HistGradientBoostingClassifier

FEAT = "/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED  = sys.argv[1] if len(sys.argv) > 1 else "Y130G-clone2"
CAL = sys.argv[2] if len(sys.argv) > 2 else "nCas9-clone1"
MASK, SEED = "Parent", 0
FEATS = ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

d = np.load(f"{FEAT}/a3a_trainset_v5.npz", allow_pickle=True)
Xtr = np.column_stack([d[f] for f in FEATS]).astype(np.float32)
ytr = d["y"].astype(np.int8)
log(f"PCAWG train {Xtr.shape}, {int(ytr.sum()):,} positives, base rate {ytr.mean():.4f}")
m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True,
                                   validation_fraction=0.1, random_state=SEED).fit(Xtr, ytr)
log("model fitted on PCAWG ONLY -- no HEK293T data has been seen")

sc, ed_s, cal_s = [], [], []
for c in CH:
    U = np.load(f"{FEAT}/universe_chr{c}.npz")
    E = np.load(f"{FEAT}/counts_{ED}_chr{c}.npz")
    K = np.load(f"{FEAT}/counts_{CAL}_chr{c}.npz")
    P = np.load(f"{FEAT}/counts_{MASK}_chr{c}.npz")
    el = (E["cov"] >= 8) & (K["cov"] >= 8) & (P["cov"] >= 15) & (P["alt"] == 0)
    X = np.column_stack([U[f][el] for f in FEATS]).astype(np.float32)
    sc.append(m.predict_proba(X)[:, 1].astype(np.float32))
    ed_s.append(((E["alt"] >= 2) & (K["alt"] == 0))[el])
    cal_s.append(((K["alt"] >= 2) & (E["alt"] == 0))[el])
    del U, E, K, P, X
sc = np.concatenate(sc); ed = np.concatenate(ed_s); cal = np.concatenate(cal_s)
N = len(sc)
log(f"scored {N:,} jointly-eligible sites   editor-specific {int(ed.sum()):,}   "
    f"calibrator-specific {int(cal.sum()):,}")

order = np.argsort(-sc)
rng = np.random.default_rng(SEED)
print(f"\n  {'panel':>8} {'n_sites':>12} {'ED captured':>12} {'cap%':>7} {'x':>7} "
      f"{'CAL captured':>13} {'cap%':>7} {'x':>7} {'ed-cal':>8}")
for k in (0.01, 0.1, 1.0, 5.0):
    n = max(1, int(N * k / 100))
    top = order[:n]
    e_cap = int(ed[top].sum()); c_cap = int(cal[top].sum())
    e_pct = 100 * e_cap / max(int(ed.sum()), 1); c_pct = 100 * c_cap / max(int(cal.sum()), 1)
    e_x = e_pct / k; c_x = c_pct / k
    print(f"  {k:7.2f}% {n:12,} {e_cap:12,} {e_pct:6.2f}% {e_x:7.3f} "
          f"{c_cap:13,} {c_pct:6.2f}% {c_x:7.3f} {e_x-c_x:+8.3f}")
# RANDOM baseline beside every enrichment, as required
print("\n  RANDOM baseline (same panel sizes, shuffled ranking):")
for k in (0.1, 1.0, 5.0):
    n = max(1, int(N * k / 100))
    r = rng.permutation(N)[:n]
    e_x = (100 * int(ed[r].sum()) / max(int(ed.sum()), 1)) / k
    c_x = (100 * int(cal[r].sum()) / max(int(cal.sum()), 1)) / k
    print(f"  {k:7.2f}%   editor {e_x:6.3f}x   calibrator {c_x:6.3f}x")
print("\n  A panel capturing K% of mutations from K% of sites is 1.000x = NO SKILL.")
print("  Only an editor-minus-calibrator gap is transfer.")
