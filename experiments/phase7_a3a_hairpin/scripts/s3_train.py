#!/usr/bin/env python
"""
S3 — A3A endogenous classifier, trained on PCAWG, held out by CHROMOSOME.

The scientific question is NOT "can we get a good AUROC" -- with trinuc-matched
negatives even a mediocre AUROC is meaningful, and in this project a high AUROC
has repeatedly meant motif-learning rather than site prediction. The questions are:

  1. Does ANY beyond-motif signal exist?          (AUROC vs 0.5, held-out chrom)
  2. Does explicit STRUCTURE add over sequence?   (ablation: hp / seq / both)
  3. Does it CONCENTRATE?                         (enrichment at top-K%, with
                                                   RANDOM baseline printed)

Ablation is the point: extended sequence context can partially encode
hairpin-ness on its own. If hp-only ~ seq-only ~ both, the structure features are
redundant re-descriptions of local sequence, not an independent mechanism.
"""
import os, time, json
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

BASE = "/mnt/data/a3a"
FEAT = f"{BASE}/feat"
REF = "/mnt/data/ref/hg19.fa"
CTX = 10          # +/- bp of extended sequence context
SEED = 20260810
rng = np.random.default_rng(SEED)
B = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_hg19():
    seqs, name, chunks = {}, None, []
    with open(REF) as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0].replace("chr", "")
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks).upper()
    out = {}
    for c in [str(i) for i in range(1, 23)] + ["X"]:
        if c in seqs:
            a = np.frombuffer(seqs[c].encode(), dtype=np.uint8)
            o = np.full(a.shape, 4, dtype=np.uint8)
            for ch, v in B.items():
                if ch != "N":
                    o[a == ord(ch)] = v
            out[c] = o
    return out


def seq_context(genome, chrom, pos0):
    """One-hot +/-CTX bp, focal base excluded (always C)."""
    n = len(pos0)
    offs = np.array([o for o in range(-CTX, CTX + 1) if o != 0])
    codes = np.full((n, len(offs)), 4, dtype=np.uint8)
    for c in np.unique(chrom):
        m = chrom == c
        seq = genome[c]
        idx = np.clip(pos0[m][:, None] + offs[None, :], 0, len(seq) - 1)
        codes[m] = seq[idx]
    oh = np.zeros((n, len(offs) * 4), dtype=np.float32)
    for j in range(len(offs)):
        for b in range(4):
            oh[:, j * 4 + b] = (codes[:, j] == b)
    names = [f"seq{o:+d}_{'ACGT'[b]}" for o in offs for b in range(4)]
    return oh, names


def enrichment(score, y, label):
    """Top-K% enrichment of positives, with a shuffled-score RANDOM baseline."""
    order = np.argsort(-score)
    base = y.mean()
    out = []
    for k in [0.1, 1.0, 5.0, 10.0]:
        n = max(int(len(y) * k / 100), 1)
        e = y[order[:n]].mean() / base
        sh = rng.permutation(score)
        r = y[np.argsort(-sh)[:n]].mean() / base
        out.append((k, e, r))
        log(f"    {label} top{k:>5.1f}%  enr={e:6.3f}x   RANDOM={r:5.3f}x   (n={n:,})")
    return out


def main():
    t0 = time.time()
    d = np.load(f"{FEAT}/a3a_trainset.npz", allow_pickle=True)
    y = d["y"].astype(np.int8)
    chrom = d["chrom"]
    pos0 = d["pos"].astype(np.int64) - 1
    log(f"trainset: {len(y):,} sites, {int(y.sum()):,} positive")

    hp = np.column_stack([d["stem"], d["loop"], d["kpos"], d["gc_pairs"],
                          d["hp_score"], d["local_gc"]]).astype(np.float32)
    hp_names = ["stem", "loop", "pos_in_loop", "gc_pairs", "hp_score", "local_gc"]

    log("building extended sequence context")
    genome = load_hg19()
    sq, sq_names = seq_context(genome, chrom, pos0)
    del genome

    blocks = {"hairpin_only": (hp, hp_names),
              "sequence_only": (sq, sq_names),
              "hairpin+sequence": (np.hstack([hp, sq]), hp_names + sq_names)}

    # held-out CHROMOSOME folds (never random -- spatial leakage inflated a
    # previous result in this project)
    chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
    folds = [chroms[i::5] for i in range(5)]
    log(f"5 folds by chromosome: {[len(f) for f in folds]} chroms each")

    results = {}
    for bname, (X, names) in blocks.items():
        log(f"=== {bname}  ({X.shape[1]} features)")
        aucs, oof = [], np.zeros(len(y))
        for fi, hold in enumerate(folds):
            te = np.isin(chrom, hold)
            tr = ~te
            clf = HistGradientBoostingClassifier(
                max_iter=300, learning_rate=0.08, max_leaf_nodes=63,
                min_samples_leaf=100, l2_regularization=1.0,
                early_stopping=True, validation_fraction=0.1,
                random_state=SEED)
            clf.fit(X[tr], y[tr])
            p = clf.predict_proba(X[te])[:, 1]
            oof[te] = p
            a = roc_auc_score(y[te], p)
            aucs.append(a)
            log(f"  fold{fi} heldout={len(hold)}chr n_te={te.sum():,} AUROC={a:.4f}")
        log(f"  {bname}: AUROC mean={np.mean(aucs):.4f} sd={np.std(aucs):.4f}")
        enr = enrichment(oof, y, bname)
        results[bname] = dict(auc_mean=float(np.mean(aucs)),
                              auc_sd=float(np.std(aucs)),
                              auc_folds=[float(a) for a in aucs],
                              enrichment=[[float(x) for x in e] for e in enr])
        np.save(f"{FEAT}/oof_{bname}.npy", oof)

    log("=" * 64)
    log("ABLATION SUMMARY (held-out chromosome, trinuc-matched negatives)")
    for k, v in results.items():
        log(f"  {k:20s} AUROC={v['auc_mean']:.4f}+-{v['auc_sd']:.4f}  "
            f"top1%={v['enrichment'][1][1]:.3f}x (rand {v['enrichment'][1][2]:.3f}x)")
    hs = results["hairpin+sequence"]["auc_mean"]
    so = results["sequence_only"]["auc_mean"]
    ho = results["hairpin_only"]["auc_mean"]
    log(f"  structure adds over sequence: {hs - so:+.4f} AUROC")
    log(f"  sequence adds over structure: {hs - ho:+.4f} AUROC")
    log("  NOTE negatives trinuc-matched => none of this is motif-learning.")
    log("=" * 64)
    json.dump(results, open(f"{FEAT}/s3_results.json", "w"), indent=2)
    open(f"{BASE}/flags/S3_TRAIN_DONE", "w").write("ok\n")
    log(f"S3 complete in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
