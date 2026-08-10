#!/usr/bin/env python
"""
S5 — context-length sweep. Cheap CPU proxy for "would a DNA LM help?"

The ablation showed a +/-10bp sequence model reaches AUROC 0.618 and explicit
hairpin features add ~nothing on top (+0.0024). The natural next question is
whether A3A targeting is sequence-encodable BEYOND +/-10bp -- i.e. the
information-floor question this program keeps running into.

A DNA LM is the heavyweight way to ask that. This is the cheap way: sweep the
context window in the SAME gradient-boosting model and see whether performance is
still climbing at +/-10bp or has already plateaued.

  still climbing at the top of the sweep -> longer-range signal exists, a DNA LM
                                            is worth the setup cost
  plateaued (or worse) by ~+/-20bp        -> +/-10bp is already the ceiling and an
                                            LM will not rescue it; do not spend
                                            the GPU on it

Uses the strand-oriented windows stored by stage2b, so no genome re-read and no
opportunity to reintroduce an orientation bug. Held out by CHROMOSOME.
"""
import time, json
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

FEAT = "/mnt/data/a3a/feat"
CTXS = [1, 2, 3, 5, 10, 15, 20, 30, 40]
SEED = 20260810
rng = np.random.default_rng(SEED)


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main():
    d = np.load(f"{FEAT}/a3a_trainset_v2.npz", allow_pickle=True)
    y = d["y"].astype(np.int8); chrom = d["chrom"]; win = d["win"]
    F = win.shape[1] // 2
    log(f"{len(y):,} sites, {int(y.sum()):,} positive, window half-width {F}")

    chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
    folds = [chroms[i::5] for i in range(5)]
    out = {}
    for ctx in CTXS:
        if ctx > F:
            continue
        offs = [o for o in range(-ctx, ctx + 1) if o != 0]
        X = np.zeros((len(y), len(offs) * 4), dtype=np.float32)
        for j, o in enumerate(offs):
            col = win[:, F + o]
            for b in range(4):
                X[:, j * 4 + b] = (col == b)
        aucs, oof = [], np.zeros(len(y))
        for hold in folds:
            te = np.isin(chrom, hold); tr = ~te
            clf = HistGradientBoostingClassifier(
                max_iter=300, learning_rate=0.08, max_leaf_nodes=63,
                min_samples_leaf=100, l2_regularization=1.0,
                early_stopping=True, validation_fraction=0.1, random_state=SEED)
            clf.fit(X[tr], y[tr])
            p = clf.predict_proba(X[te])[:, 1]
            oof[te] = p
            aucs.append(roc_auc_score(y[te], p))
        n = max(int(len(y) * 0.01), 1)
        top = np.argsort(-oof)[:n]
        enr = y[top].mean() / y.mean()
        sh = rng.permutation(oof)
        rnd = y[np.argsort(-sh)[:n]].mean() / y.mean()
        out[ctx] = dict(auc=float(np.mean(aucs)), sd=float(np.std(aucs)),
                        top1=float(enr), rand=float(rnd), nfeat=X.shape[1])
        log(f"  ctx +/-{ctx:>2}bp ({X.shape[1]:>3} feat): AUROC={np.mean(aucs):.4f}"
            f"+-{np.std(aucs):.4f}  top1%={enr:6.3f}x  RANDOM={rnd:5.3f}x")
        del X

    log("=" * 62)
    ks = sorted(out)
    log("CONTEXT SWEEP (held-out chromosome, strand+trinuc-matched negatives)")
    for i, k in enumerate(ks):
        delta = "" if i == 0 else f"  d={out[k]['auc']-out[ks[i-1]]['auc']:+.4f}"
        log(f"  +/-{k:>2}bp  AUROC={out[k]['auc']:.4f}  top1%={out[k]['top1']:.3f}x{delta}")
    best = max(ks, key=lambda k: out[k]["auc"])
    log(f"  best context = +/-{best}bp (AUROC {out[best]['auc']:.4f})")
    tail = [out[k]["auc"] for k in ks if k >= 10]
    if len(tail) >= 2 and max(tail) - tail[0] < 0.005:
        log("  VERDICT: plateaued by +/-10bp => a DNA LM is NOT worth the GPU setup;")
        log("           local sequence is already saturated (information floor).")
    else:
        log("  VERDICT: still climbing => longer-range signal exists; a DNA LM IS")
        log("           worth the setup cost.")
    log("=" * 62)
    json.dump(out, open(f"{FEAT}/s5_ctxsweep.json", "w"), indent=2)
    open("/mnt/data/a3a/flags/S5_CTXSWEEP_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
