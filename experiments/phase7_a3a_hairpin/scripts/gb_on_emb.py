#!/usr/bin/env python
"""Encoder or head? The comparison that makes the DeaminaFormer verdict interpretable.

The flat-MLP ladder scored 2.631x where GB scored 5.280x on the same folds and metric. Two
readings survive that:
    (a) the ENCODER matters -- one-hot needs a convolution, and a foundation model should
        beat both;
    (b) the HEAD matters -- gradient boosting is simply a better learner on this problem, at
        any input.
The CNN ladder tests (a). This tests (b), and without it a foundation-model win or loss is
uninterpretable: if the embeddings beat GB-on-one-hot under an MLP head we would credit the
encoder, when the head might be doing the work -- or hiding it.

So: the SAME gradient booster that produced 5.280x, run on the foundation-model embeddings
and the structure blocks. Reduced by PCA fitted on the TRAINING FOLD ONLY -- fitting PCA on
all rows would leak held-out chromosomes into the projection, which is the quiet version of
the leak this project keeps finding.

Same folds, same pooled-OOF top-K% metric, random baseline beside every number.
"""
import os, json, time
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.decomposition import PCA

BASE = "/mnt/a3a"
NPC = int(os.environ.get("GB_PCA", "128"))
tr = np.load(f"{BASE}/a3a_trainset_v5.npz", allow_pickle=True)
y = tr["y"].astype(np.int8); chrom = tr["chrom"]; N = len(y)
base = float(y.mean())
print(f"v5 n={N:,} base rate {base:.6f} ceiling {1/base:.3f}x", flush=True)

STRUCT = np.stack([tr[k].astype(np.float32) for k in
                   ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")], 1)

def dense(p, key=None):
    if not os.path.exists(p):
        print(f"  {os.path.basename(p)} absent -- skipped", flush=True); return None
    z = np.load(p, allow_pickle=True)
    if len(z["y"]) != N or not (z["pos"] == tr["pos"]).all():
        print(f"  {os.path.basename(p)} row mismatch -- SKIPPED", flush=True); return None
    return (np.concatenate([z["mean"], z["cls"]], 1).astype(np.float32) if key is None
            else z[key].astype(np.float32))

NT = dense(f"{BASE}/emb_ntv2_v5.npz")
HY = dense(f"{BASE}/emb_hyena_v5.npz")
TH = dense(f"{BASE}/struct_thermo_v5.npz", "X")

chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
folds = [chroms[i::5] for i in range(5)]

def tail(score, pct):
    k = max(int(round(N * pct / 100)), 1)
    idx = np.argpartition(-score, k - 1)[:k]
    return float(y[idx].mean() / base), k

def run(name, blocks, pca_blocks):
    oof = np.zeros(N, np.float32)
    for hold in folds:
        te = np.isin(chrom, hold); trn = ~te
        parts_tr, parts_te = [], []
        for B, do_pca in blocks:
            if do_pca:
                # PCA fitted on the TRAINING FOLD ONLY -- never on all rows
                p = PCA(n_components=min(NPC, B.shape[1]), random_state=0)
                sub = np.flatnonzero(trn)
                sub = np.random.default_rng(0).choice(sub, min(150000, len(sub)), replace=False)
                p.fit(B[sub])
                parts_tr.append(p.transform(B[trn])); parts_te.append(p.transform(B[te]))
            else:
                parts_tr.append(B[trn]); parts_te.append(B[te])
        Xtr = np.concatenate(parts_tr, 1); Xte = np.concatenate(parts_te, 1)
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1,
                                           early_stopping=True, validation_fraction=0.1,
                                           random_state=0)
        m.fit(Xtr, y[trn])
        oof[te] = m.predict_proba(Xte)[:, 1]
        del Xtr, Xte
    rnd = np.random.default_rng(0).permutation(N).astype(np.float32)
    e01, k01 = tail(oof, 0.1); r01, _ = tail(rnd, 0.1)
    e1, _ = tail(oof, 1.0); r1, _ = tail(rnd, 1.0)
    print(f"  {name:34s} top0.1% {e01:6.3f}x (rand {r01:.3f}, n={k01:,})  "
          f"top1% {e1:6.3f}x (rand {r1:.3f})", flush=True)
    return {"top0.1": e01, "rand0.1": r01, "n": k01, "top1.0": e1}

L = [("GB struct only", [(STRUCT, False)])]
if TH is not None: L.append(("GB struct+thermo", [(STRUCT, False), (TH, False)]))
if NT is not None:
    L.append((f"GB ntv2(PCA{NPC})", [(NT, True)]))
    L.append((f"GB ntv2(PCA{NPC})+struct", [(NT, True), (STRUCT, False)]))
if NT is not None and HY is not None and TH is not None:
    L.append(("GB ALL(PCA)+struct+thermo",
              [(NT, True), (HY, True), (STRUCT, False), (TH, False)]))

res = {}
for name, blocks in L:
    t0 = time.time()
    res[name] = run(name, blocks, None)
    res[name]["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(res, open(f"{BASE}/gb_on_embeddings.json", "w"), indent=2)
print("\n  GB on one-hot+hairpin (the published baseline): 5.280x")
print("  flat-MLP on the same one-hot+hairpin:            2.631x")
print("  If GB-on-embeddings lands near 5.280x, the HEAD explained the MLP gap, not the")
print("  encoder. If it lands well above, the embeddings carry something new.")
