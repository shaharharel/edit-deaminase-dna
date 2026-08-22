#!/usr/bin/env python
"""DeaminaFormer-DNA, stage 2: the multi-modal fusion and its ablation ladder.

The question is not "does a foundation model help" -- it is WHICH MODALITY CARRIES THE TAIL.
The GB baseline already showed the answer is not more sequence: 80 one-hot features gave
3.814x while six hairpin-geometry features gave 5.017x, and fusing them gave 5.280x. So the
ladder below is built to locate the information, not to produce one headline number.

Every block is trained identically -- same folds, same head, same schedule -- so a difference
between blocks is a difference in INPUT, which is the only thing being compared.

SELECTION RULE, fixed before any result is seen: blocks are ranked by TOP-0.1% TAIL
ENRICHMENT on held-out chromosomes. AUROC is printed because it dissociates from the tail
(the GB run had the worst-AUROC block nearly winning the tail) and printing it keeps that
visible -- it is never the selection criterion.

RANDOM BASELINE beside every number. Base rate 1/11, so 11.000x is the arithmetic ceiling;
anything at the ceiling is leakage, not skill.
"""
import os, sys, time, json
import numpy as np
import torch
import torch.nn as nn

BASE = "/mnt/a3a"
DEV = "cuda"
EPOCHS = int(os.environ.get("FUSE_EPOCHS", "12"))
SEED = 0
torch.manual_seed(SEED); np.random.seed(SEED)

tr = np.load(f"{BASE}/a3a_trainset_v5.npz", allow_pickle=True)
y = tr["y"].astype(np.float32); chrom = tr["chrom"]; N = len(y)
base_rate = float(y.mean())
print(f"v5: n={N:,}  base rate {base_rate:.6f}  arithmetic ceiling {1/base_rate:.3f}x", flush=True)

# --- structure / engineered block: the six hairpin-geometry features plus local GC
STRUCT = np.stack([tr[k].astype(np.float32) for k in
                   ("stem", "loop", "kpos", "gc_pairs", "hp_score", "local_gc")], 1)
# --- the 81 bp one-hot window the GB model saw, flattened.
# win uses 0-3 for ACGT and 4 for N. An earlier version of this line said np.clip(w81, 0, 3),
# which SILENTLY RECODES EVERY N AS T. Only 106 bases in 74.8 M are affected, so it would
# never have shown up in a metric -- which is precisely why it belongs in this bug family
# rather than outside it. A 5-row identity maps N to an all-zero vector: absent, not guessed.
w81 = tr["win"]
_n_amb = int((w81 > 3).sum())
EYE5 = np.zeros((5, 4), dtype=np.float32); EYE5[:4] = np.eye(4, dtype=np.float32)
SEQ81 = EYE5[w81].reshape(N, -1)
print(f"  seq81: {_n_amb} ambiguous bases encoded as all-zero (never coerced to a base)",
      flush=True)
print(f"blocks on disk: struct {STRUCT.shape}  seq81 {SEQ81.shape}", flush=True)

# thermodynamic structure block (p_unpaired from a DNA partition function, two scales).
# Optional: if it has not finished, its blocks are skipped rather than the run being blocked.
THERMO = None
_tp = f"{BASE}/struct_thermo_v5.npz"
if os.path.exists(_tp):
    _t = np.load(_tp, allow_pickle=True)
    if len(_t["y"]) == N and (_t["pos"] == tr["pos"]).all():
        THERMO = _t["X"].astype(np.float32)
        print(f"  thermo: {THERMO.shape} {list(_t['names'])}", flush=True)
    else:
        print("  thermo present but row order/length differs -- SKIPPED", flush=True)
else:
    print("  thermo block not built yet -- its ablations will be skipped", flush=True)

EMB = {}
for tag in ("ntv2", "hyena"):
    p = f"{BASE}/emb_{tag}_v5.npz"
    if not os.path.exists(p):
        print(f"  MISSING {p} -- {tag} blocks will be skipped", flush=True); continue
    e = np.load(p)
    if len(e["y"]) != N:
        print(f"  {p} has n={len(e['y']):,} but v5 has {N:,} -- SKIPPING, not truncating "
              f"(a silent length mismatch is how alignment bugs enter)", flush=True)
        continue
    # provenance: the embedding rows must be the SAME rows, in the same order
    assert (e["pos"] == tr["pos"]).all() and (e["chrom"] == chrom).all(), \
        f"{tag}: row order differs from v5 -- refusing to fuse misaligned features"
    EMB[tag] = np.concatenate([e["mean"].astype(np.float32), e["cls"].astype(np.float32)], 1)
    print(f"  {tag}: {EMB[tag].shape} (mean+cls, row order verified against v5)", flush=True)

BLOCKS = {
    "struct_only":              [STRUCT],
    "seq81_only":               [SEQ81],
    "seq81+struct (GB parity)": [SEQ81, STRUCT],
}
if THERMO is not None:
    BLOCKS["thermo_only"] = [THERMO]
    BLOCKS["struct+thermo"] = [STRUCT, THERMO]
if "ntv2" in EMB:
    BLOCKS["ntv2_only"] = [EMB["ntv2"]]
    BLOCKS["ntv2+struct"] = [EMB["ntv2"], STRUCT]
if "hyena" in EMB:
    BLOCKS["hyena_only"] = [EMB["hyena"]]
    BLOCKS["hyena+struct"] = [EMB["hyena"], STRUCT]
if "ntv2" in EMB and "hyena" in EMB:
    BLOCKS["ntv2+hyena"] = [EMB["ntv2"], EMB["hyena"]]
    BLOCKS["ntv2+hyena+struct"] = [EMB["ntv2"], EMB["hyena"], STRUCT]
    BLOCKS["ntv2+hyena+struct"] = [EMB["ntv2"], EMB["hyena"], STRUCT]
    BLOCKS["ALL (+seq81)"] = [EMB["ntv2"], EMB["hyena"], STRUCT, SEQ81]
if THERMO is not None and "ntv2" in EMB and "hyena" in EMB:
    BLOCKS["ntv2+hyena+struct+thermo"] = [EMB["ntv2"], EMB["hyena"], STRUCT, THERMO]
    BLOCKS["FULL (+seq81+thermo)"] = [EMB["ntv2"], EMB["hyena"], STRUCT, THERMO, SEQ81]

chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
folds = [chroms[i::5] for i in range(5)]
print(f"5 folds by held-out CHROMOSOME: {[len(f) for f in folds]}\n", flush=True)


class Head(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d, 512), nn.GELU(), nn.Dropout(0.2),
            nn.Linear(512, 128), nn.GELU(), nn.Dropout(0.1),
            nn.Linear(128, 1))
    def forward(self, x):
        return self.net(x).squeeze(-1)


def tail_enrichment(score, yy, pct):
    k = max(int(round(len(yy) * pct / 100)), 1)
    idx = np.argpartition(-score, k - 1)[:k]
    return float(yy[idx].mean() / base_rate), k


def run(name, mats):
    X = np.concatenate(mats, 1)
    oof = np.zeros(len(y), dtype=np.float32)
    aucs = []
    for f, hold in enumerate(folds):
        te = np.isin(chrom, hold); trn = ~te
        mu = X[trn].mean(0); sd = X[trn].std(0) + 1e-6
        Xtr = torch.from_numpy((X[trn] - mu) / sd).to(DEV)
        Xte = torch.from_numpy((X[te] - mu) / sd).to(DEV)
        ytr = torch.from_numpy(y[trn]).to(DEV)
        m = Head(X.shape[1]).to(DEV)
        opt = torch.optim.AdamW(m.parameters(), lr=1e-3, weight_decay=1e-4)
        pw = torch.tensor((1 - base_rate) / base_rate, device=DEV)
        lossf = nn.BCEWithLogitsLoss(pos_weight=pw)
        n = len(ytr); bs = 8192
        for ep in range(EPOCHS):
            perm = torch.randperm(n, device=DEV)
            m.train()
            for i in range(0, n, bs):
                b = perm[i:i + bs]
                opt.zero_grad()
                loss = lossf(m(Xtr[b]), ytr[b])
                loss.backward(); opt.step()
        m.eval()
        with torch.no_grad():
            s = torch.cat([m(Xte[i:i + 65536]) for i in range(0, len(Xte), 65536)])
        oof[te] = s.float().cpu().numpy()
        yt = y[te]
        o = np.argsort(-oof[te]); r = np.empty(len(yt)); r[o] = np.arange(len(yt))
        n1 = yt.sum(); n0 = len(yt) - n1
        aucs.append(float(((len(yt) - 1 - r[yt == 1]).sum() - n1 * (n1 - 1) / 2) / (n1 * n0)))
        del Xtr, Xte, ytr, m
        torch.cuda.empty_cache()
    rnd = np.random.default_rng(SEED).permutation(len(y)).astype(np.float32)
    res = {"auroc_mean": float(np.mean(aucs)), "auroc_sd": float(np.std(aucs)), "n_feat": X.shape[1]}
    line = f"  {name:26s} AUROC {np.mean(aucs):.4f}+-{np.std(aucs):.4f}  d={X.shape[1]:<5d}"
    for pct in (0.1, 1.0):
        e, k = tail_enrichment(oof, y, pct)
        r, _ = tail_enrichment(rnd, y, pct)
        res[f"top{pct}"] = e; res[f"rand{pct}"] = r; res[f"n_top{pct}"] = k
        line += f" | top{pct}% {e:6.3f}x (rand {r:.3f}, n={k:,})"
    print(line, flush=True)
    return res


out = {}
for name, mats in BLOCKS.items():
    t0 = time.time()
    out[name] = run(name, mats)
    out[name]["minutes"] = round((time.time() - t0) / 60, 1)
json.dump(out, open(f"{BASE}/deaminaformer_ablation.json", "w"), indent=2)

print("\n=== RANKED BY TOP-0.1% TAIL (the pre-registered selection rule) ===", flush=True)
for k, v in sorted(out.items(), key=lambda kv: -kv[1]["top0.1"]):
    print(f"  {v['top0.1']:6.3f}x  (rand {v['rand0.1']:.3f})  AUROC {v['auroc_mean']:.4f}  {k}")
print(f"\n  GB baseline to beat: 5.280x (v5) / 5.096x (coverage-matched), ceiling "
      f"{1/base_rate:.3f}x")
