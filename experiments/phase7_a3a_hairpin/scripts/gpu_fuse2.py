#!/usr/bin/env python
"""DeaminaFormer-DNA: the architecture the first ablation actually argued for.

WHAT THE FIRST RUN SHOWED. Concatenating a flattened one-hot window into an MLP is a bad
encoder for DNA: on identical folds and metric, GB scored 3.814x on the 81 bp one-hot and
the MLP scored 2.012x; fused, GB 5.280x vs MLP 2.631x. Trees split on individual positions;
a dense first layer smears them. AUROC even ROSE while the tail collapsed.

The lesson is not "MLPs are bad" -- it is that CONCATENATION IS ONLY VALID FOR DENSE
REPRESENTATIONS. So every modality gets an encoder appropriate to its own structure, and
only their dense outputs are concatenated:

    sequence, local     Conv1d motif scanner over +-200 bp   (DeepBind/Basset inductive
                                                              bias: position-invariant
                                                              motifs, which is what a
                                                              deaminase target site IS)
    sequence, global    frozen NT-v2 / HyenaDNA embeddings   (1 kb context, genome-scale
                                                              pretraining)
    structure           hairpin geometry + DNA partition     (p_unpaired is the mechanism:
                        function                              A3A needs a single-stranded C)

+-200 bp is deliberate: 2.5x GB's window, so a CNN win cannot be dismissed as "it just saw
more sequence", while staying far under the 1 kb the foundation models already cover.

SELECTION RULE, unchanged and fixed in advance: TOP-0.1% TAIL ENRICHMENT on held-out
chromosomes, against GB's 5.280x. AUROC is printed only because it dissociates. Random
baseline beside every number; base rate 1/11 makes 11.000x the arithmetic ceiling.
"""
import os, sys, time, json
import numpy as np
import torch
import torch.nn as nn

BASE = "/mnt/a3a"
DEV = "cuda"
HALF = int(os.environ.get("CNN_HALF", "200"))
EPOCHS = int(os.environ.get("FUSE_EPOCHS", "8"))
BS = int(os.environ.get("FUSE_BS", "1024"))
torch.manual_seed(0); np.random.seed(0)

# SMOKE: FUSE_LIMIT subsamples rows so a runtime error surfaces in a minute instead of at
# 3am, four hours into the ladder. Never used for a reported number -- the row count is
# printed with every result so a subsampled run cannot be mistaken for a real one.
LIMIT = int(os.environ.get("FUSE_LIMIT", "0"))
tr = dict(np.load(f"{BASE}/a3a_trainset_v5.npz", allow_pickle=True))
if LIMIT:
    keep = np.random.default_rng(0).choice(len(tr["y"]), LIMIT, replace=False); keep.sort()
    tr = {k: v[keep] for k, v in tr.items()}
y = tr["y"].astype(np.float32); chrom = tr["chrom"]; N = len(y)
base_rate = float(y.mean())
print(f"v5: n={N:,} base rate {base_rate:.6f} ceiling {1/base_rate:.3f}x", flush=True)

W = np.load(f"{BASE}/a3a_trainset_v5_win1k.npz", allow_pickle=True)
_w = W["win1k"]; _wp = W["pos"]
if LIMIT:
    _w = _w[keep]; _wp = _wp[keep]
assert (_wp == tr["pos"]).all(), "win1k row order differs from v5"
MID = _w.shape[1] // 2
SEQ = np.ascontiguousarray(_w[:, MID - HALF:MID + HALF + 1])   # uint8 0-4
print(f"CNN input: {SEQ.shape} (+-{HALF} bp, vs GB's +-40)", flush=True)

STRUCT = np.stack([tr[k].astype(np.float32) for k in
                   ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")], 1)

def load_dense(path, key=None):
    if not os.path.exists(path):
        print(f"  {os.path.basename(path)} absent -- its blocks are skipped", flush=True)
        return None
    z = dict(np.load(path, allow_pickle=True))
    if LIMIT and len(z["y"]) == len(keep) + 0 or (LIMIT and len(z["y"]) > N):
        z = {k: (v[keep] if hasattr(v, "__len__") and len(v) > N else v) for k, v in z.items()}
    if len(z["y"]) != N or not (z["pos"] == tr["pos"]).all():
        print(f"  {os.path.basename(path)} row order/length mismatch -- SKIPPED, not "
              f"truncated", flush=True)
        return None
    return (np.concatenate([z["mean"], z["cls"]], 1).astype(np.float32) if key is None
            else z[key].astype(np.float32))

NT = load_dense(f"{BASE}/emb_ntv2_v5.npz")
HY = load_dense(f"{BASE}/emb_hyena_v5.npz")
TH = load_dense(f"{BASE}/struct_thermo_v5.npz", "X")

DENSE = {"struct": STRUCT}
if TH is not None: DENSE["thermo"] = TH
if NT is not None: DENSE["ntv2"] = NT
if HY is not None: DENSE["hyena"] = HY
print(f"  dense blocks available: {list(DENSE)}", flush=True)


class MotifCNN(nn.Module):
    """Position-invariant motif detectors, the inductive bias a flattened MLP lacks."""
    def __init__(self, out=128):
        super().__init__()
        self.body = nn.Sequential(
            nn.Conv1d(5, 128, 11, padding=5), nn.ReLU(), nn.MaxPool1d(4),
            nn.Conv1d(128, 128, 7, padding=3), nn.ReLU(), nn.MaxPool1d(4),
            nn.Conv1d(128, 64, 5, padding=2), nn.ReLU())
        self.proj = nn.Linear(128, out)
    def forward(self, x):                       # x int64 [B, L]
        h = self.body(torch.nn.functional.one_hot(x, 5).permute(0, 2, 1).float())
        h = torch.cat([h.max(-1).values, h.mean(-1)], 1)
        return self.proj(h)


class Net(nn.Module):
    def __init__(self, use_cnn, dense_dim):
        super().__init__()
        self.cnn = MotifCNN() if use_cnn else None
        d = (128 if use_cnn else 0) + (128 if dense_dim else 0)
        self.dense = nn.Sequential(nn.Linear(dense_dim, 256), nn.GELU(),
                                   nn.Linear(256, 128), nn.GELU()) if dense_dim else None
        self.head = nn.Sequential(nn.Linear(d, 128), nn.GELU(), nn.Dropout(0.15),
                                  nn.Linear(128, 1))
    def forward(self, seq, dz):
        parts = []
        if self.cnn is not None: parts.append(self.cnn(seq))
        if self.dense is not None: parts.append(self.dense(dz))
        return self.head(torch.cat(parts, 1)).squeeze(-1)


chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
folds = [chroms[i::5] for i in range(5)]
print(f"5 folds by held-out CHROMOSOME: {[len(f) for f in folds]}\n", flush=True)


def tail(score, pct):
    k = max(int(round(N * pct / 100)), 1)
    idx = np.argpartition(-score, k - 1)[:k]
    return float(y[idx].mean() / base_rate), k


def run(name, use_cnn, dense_keys):
    Xd = np.concatenate([DENSE[k] for k in dense_keys], 1) if dense_keys else None
    oof = np.zeros(N, np.float32); aucs = []
    # keep the window as uint8 (371 MB) and cast per batch on the GPU. Casting the
    # whole array to int64 up front costs 2.96 GB and was being redone for every
    # block in the ladder.
    seq_t = torch.from_numpy(SEQ)
    for hold in folds:
        te = np.isin(chrom, hold); trn = ~te
        if Xd is not None:
            mu, sd = Xd[trn].mean(0), Xd[trn].std(0) + 1e-6
        m = Net(use_cnn, Xd.shape[1] if Xd is not None else 0).to(DEV)
        opt = torch.optim.AdamW(m.parameters(), lr=1e-3, weight_decay=1e-4)
        lossf = nn.BCEWithLogitsLoss(
            pos_weight=torch.tensor((1 - base_rate) / base_rate, device=DEV))
        idx_tr = np.flatnonzero(trn); idx_te = np.flatnonzero(te)
        for ep in range(EPOCHS):
            m.train()
            perm = np.random.permutation(idx_tr)
            for i in range(0, len(perm), BS):
                b = perm[i:i + BS]
                s = seq_t[b].to(DEV, non_blocking=True).long()
                dz = (torch.from_numpy((Xd[b] - mu) / sd).to(DEV)
                      if Xd is not None else None)
                opt.zero_grad()
                lossf(m(s, dz), torch.from_numpy(y[b]).to(DEV)).backward()
                opt.step()
        m.eval(); out = []
        with torch.no_grad():
            for i in range(0, len(idx_te), 4096):
                b = idx_te[i:i + 4096]
                s = seq_t[b].to(DEV).long()
                dz = (torch.from_numpy((Xd[b] - mu) / sd).to(DEV)
                      if Xd is not None else None)
                out.append(m(s, dz).float().cpu().numpy())
        oof[te] = np.concatenate(out)
        yt = y[te]; o = np.argsort(-oof[te]); r = np.empty(len(yt)); r[o] = np.arange(len(yt))
        n1 = yt.sum(); n0 = len(yt) - n1
        aucs.append(float(((len(yt) - 1 - r[yt == 1]).sum() - n1 * (n1 - 1) / 2) / (n1 * n0)))
        del m; torch.cuda.empty_cache()
    rnd = np.random.default_rng(0).permutation(N).astype(np.float32)
    e01, k01 = tail(oof, 0.1); r01, _ = tail(rnd, 0.1)
    e1, k1 = tail(oof, 1.0); r1, _ = tail(rnd, 1.0)
    print(f"  [n={N:,}] {name:30s} AUROC {np.mean(aucs):.4f}  top0.1% {e01:6.3f}x (rand {r01:.3f}, "
          f"n={k01:,})  top1% {e1:6.3f}x (rand {r1:.3f})", flush=True)
    np.save(f"{BASE}/oof_fuse2_{name.replace(' ', '_')}.npy", oof)
    del Xd, seq_t
    import gc; gc.collect()
    return {"auroc": float(np.mean(aucs)), "top0.1": e01, "rand0.1": r01, "n_top": k01,
            "top1.0": e1}


LADDER = [("cnn_only", True, []), ("struct_only(MLP)", False, ["struct"]),
          ("cnn+struct", True, ["struct"])]
if "thermo" in DENSE:
    LADDER += [("thermo_only(MLP)", False, ["thermo"]),
               ("cnn+struct+thermo", True, ["struct", "thermo"])]
if "ntv2" in DENSE:
    LADDER += [("ntv2_only(MLP)", False, ["ntv2"]), ("cnn+ntv2", True, ["ntv2"])]
if "hyena" in DENSE:
    LADDER += [("hyena_only(MLP)", False, ["hyena"])]
ALL = [k for k in ("ntv2", "hyena", "struct", "thermo") if k in DENSE]
if len(ALL) >= 3:
    LADDER.append(("DEAMINAFORMER (cnn+all)", True, ALL))

res = {}
for name, uc, dk in LADDER:
    t0 = time.time()
    res[name] = run(name, uc, dk)
    res[name]["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(res, open(f"{BASE}/deaminaformer_v2.json", "w"), indent=2)

print("\n=== RANKED BY TOP-0.1% TAIL (pre-registered rule) ===", flush=True)
for k, v in sorted(res.items(), key=lambda kv: -kv[1]["top0.1"]):
    print(f"  {v['top0.1']:6.3f}x (rand {v['rand0.1']:.3f})  AUROC {v['auroc']:.4f}  {k}")
print(f"\n  GB baseline: 5.280x   flat-MLP fusion: 2.631x   ceiling {1/base_rate:.3f}x")
