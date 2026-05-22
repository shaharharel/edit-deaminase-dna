"""DeaminaFormer-DNA training entry on 74,759-site combined parquet.

Usage:
    python experiments/phase2_train/train_dna.py \
        --embeddings ~/dna_emb/hyenadna_dna_offtarget_combined.pt \
        --emb-dim 256 \
        --encode data/processed/encode_features.parquet \
        --device cuda

Or Evo-2:
    --embeddings ~/dna_emb/evo2_dna_offtarget.pt --emb-dim 1024

LOO CV by `source_paper` so we test cross-source generalization (not random fold).
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve


def focal_loss(logits, targets, gamma=2.0, alpha=None):
    """Binary focal loss. alpha = positive class weight; if None, no rebalancing."""
    bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
    p = torch.sigmoid(logits)
    pt = p * targets + (1 - p) * (1 - targets)
    focal_term = (1 - pt) ** gamma
    if alpha is not None:
        alpha_t = alpha * targets + (1 - alpha) * (1 - targets)
        focal_term = focal_term * alpha_t
    return (focal_term * bce).mean()


def recall_at_fpr(y_true, y_score, fpr_target=0.01):
    """Recall at a target false-positive rate."""
    from sklearn.metrics import roc_curve
    fpr, tpr, _ = roc_curve(y_true, y_score)
    # Interpolate to find recall at the exact FPR target
    if fpr[-1] < fpr_target:
        return float(tpr[-1])
    return float(np.interp(fpr_target, fpr, tpr))

sys.path.insert(0, "/Users/shaharharel/Documents/github/edit-deaminase-dna")
sys.path.insert(0, "/Users/shaharharel/Documents/github/edit-deaminase")
from src.models.deaminaformer_dna import DeaminaFormerDNA, DNA_SUBSTRATES  # type: ignore
from src.data.editor_features import CATALOG  # type: ignore


def _site_key(r):
    return f"{r['chrom']}:{int(r['pos'])}:{r['strand']}:{r['ref']}:{r['alt']}"


def _editor_lookup(name: str):
    if name in CATALOG:
        return CATALOG[name]
    # Try common aliases
    if "ABE" in name.upper():
        return CATALOG.get("ABE8e", next(iter(CATALOG.values())))
    if "BE" in name.upper() or "CBE" in name.upper():
        return CATALOG.get("BE4max", next(iter(CATALOG.values())))
    return CATALOG.get("UNK", next(iter(CATALOG.values())))


class DNADataset(Dataset):
    def __init__(self, df: pd.DataFrame, embeddings: dict, encode_df: pd.DataFrame | None,
                  emb_dim: int):
        self.df = df.reset_index(drop=True)
        self.emb = embeddings
        self.emb_dim = emb_dim
        self.encode = ({(r["chrom"], int(r["pos"]), r["strand"]):
                         (r["dnase_mean"], r["rloop_mean"], r["mappability_mean"],
                          r["atac_overlap"]) for _, r in encode_df.iterrows()}
                        if encode_df is not None else {})

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        r = self.df.iloc[idx]
        key = _site_key(r)
        emb = self.emb.get(key)
        if emb is None:
            emb = torch.zeros(self.emb_dim, dtype=torch.float32)
        emb = emb.float()
        enc_key = (r["chrom"], int(r["pos"]), r["strand"])
        enc = self.encode.get(enc_key, (0.0, 0.0, 1.0, 0))
        # ENCODE features: dnase, rloop, mappability, atac_overlap + GC + region_class
        # GC placeholder (computed on the fly later); region_class = 0 for now
        encode_feats = np.array([enc[0], enc[1], enc[2], float(enc[3]), 0.5, 0.0], dtype=np.float32)
        # Substrate token: open chromatin if ATAC+, else replication fork
        sub = 2 if enc[3] else 1   # 2=open_chromatin, 1=replication_fork
        return {
            "emb": emb,
            "encode": encode_feats,
            "editor": str(r["editor"]),
            "substrate_id": sub,
            "label": float(r["is_positive"]),
        }


def collate(batch):
    emb = torch.stack([b["emb"] for b in batch])
    enc = torch.from_numpy(np.stack([b["encode"] for b in batch]))
    eds = [_editor_lookup(b["editor"]) for b in batch]
    subs = torch.tensor([b["substrate_id"] for b in batch], dtype=torch.long)
    y = torch.tensor([b["label"] for b in batch], dtype=torch.float32)
    return {"emb": emb, "enc": enc, "editors": eds, "substrate_ids": subs, "labels": y}


def train_one_fold(train_df, val_df, *, emb_dict, encode_df, emb_dim, device, epochs, lr, batch_size, patience=4):
    train_ds = DNADataset(train_df, emb_dict, encode_df, emb_dim)
    val_ds = DNADataset(val_df, emb_dict, encode_df, emb_dim)
    train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True, collate_fn=collate)
    val_dl = DataLoader(val_ds, batch_size=batch_size * 4, shuffle=False, collate_fn=collate)

    model = DeaminaFormerDNA(d_evo2=emb_dim).to(device)
    optim = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-3)

    best_auroc = 0.0
    best_auprc = 0.0
    bad = 0
    for ep in range(epochs):
        model.train()
        ep_loss = 0.0; n = 0
        for batch in train_dl:
            emb = batch["emb"].to(device)
            enc = batch["enc"].to(device)
            sub = batch["substrate_ids"].to(device)
            y = batch["labels"].to(device)
            out = model(emb, enc, batch["editors"], sub)
            loss = F.binary_cross_entropy_with_logits(out["logit"], y)
            optim.zero_grad(); loss.backward(); optim.step()
            ep_loss += loss.item() * len(y); n += len(y)
        model.eval()
        ys, ps = [], []
        with torch.no_grad():
            for batch in val_dl:
                emb = batch["emb"].to(device)
                enc = batch["enc"].to(device)
                sub = batch["substrate_ids"].to(device)
                out = model(emb, enc, batch["editors"], sub)
                ys.append(batch["labels"].numpy())
                ps.append(torch.sigmoid(out["logit"]).cpu().numpy())
        ys = np.concatenate(ys); ps = np.concatenate(ps)
        try:
            auroc = roc_auc_score(ys, ps); auprc = average_precision_score(ys, ps)
        except ValueError:
            auroc = float("nan"); auprc = float("nan")
        improved = auroc > best_auroc
        print(f"  ep{ep+1}/{epochs}  loss={ep_loss/max(n,1):.4f}  val_auroc={auroc:.4f}  auprc={auprc:.4f}{'  *' if improved else ''}")
        if improved:
            best_auroc, best_auprc, bad = auroc, auprc, 0
        else:
            bad += 1
            if bad >= patience:
                print(f"  early-stop ep{ep+1}")
                break
    return {"best_auroc": best_auroc, "best_auprc": best_auprc, "n_val": len(ys)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_combined.parquet")
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--emb-dim", type=int, required=True)
    ap.add_argument("--encode", default=None, help="Optional ENCODE features parquet")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--epochs", type=int, default=15)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--out", default="experiments/phase2_train/outputs_dna/results.json")
    args = ap.parse_args()

    print(f"[train] embeddings={args.embeddings} emb_dim={args.emb_dim}")
    emb_dict = torch.load(args.embeddings, map_location="cpu")
    print(f"[train] {len(emb_dict)} embeddings loaded")

    df = pd.read_parquet(args.combined)
    print(f"[train] {len(df)} sites (pos={df.is_positive.sum()} neg={(~df.is_positive).sum()})")

    encode_df = pd.read_parquet(args.encode) if args.encode and Path(args.encode).exists() else None
    if encode_df is not None:
        print(f"[train] ENCODE features: {len(encode_df)} sites")
    else:
        print(f"[train] (no ENCODE features — using zeros)")

    device = torch.device(args.device)

    # LOO CV: leave one positive source_paper out (+ that source's negatives based on matched_neg)
    pos_sources = sorted(df[df.is_positive].source_paper.unique())
    # Skip tiny sources (<200 positives)
    pos_sources = [s for s in pos_sources if (df[df.source_paper == s].shape[0]) >= 200]
    print(f"[train] LOO sources ({len(pos_sources)}): {pos_sources}")

    results = {}
    for holdout in pos_sources:
        val_mask = df["source_paper"] == holdout
        train_df = df[~val_mask].reset_index(drop=True)
        val_df = df[val_mask].reset_index(drop=True)
        # Need both classes in val
        if val_df.is_positive.nunique() < 2:
            # Add some matched_neg to val if val is all-positive holdout
            n_neg_needed = min(len(val_df), df[~df.is_positive].shape[0] // 4)
            extra_neg = df[~df.is_positive].sample(n=n_neg_needed, random_state=42)
            val_df = pd.concat([val_df, extra_neg], ignore_index=True)
            train_df = train_df[~train_df.index.isin(extra_neg.index)].reset_index(drop=True)
        print(f"\n=== Fold: holdout={holdout} ===")
        print(f"  train={len(train_df)} (pos {train_df.is_positive.sum()} / neg {(~train_df.is_positive).sum()})")
        print(f"  val={len(val_df)} (pos {val_df.is_positive.sum()} / neg {(~val_df.is_positive).sum()})")
        metrics = train_one_fold(train_df, val_df,
                                  emb_dict=emb_dict, encode_df=encode_df, emb_dim=args.emb_dim,
                                  device=device, epochs=args.epochs,
                                  lr=args.lr, batch_size=args.batch_size)
        results[holdout] = metrics

    out_p = Path(args.out)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    with open(out_p, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[done] Wrote {out_p}")
    for k, v in results.items():
        print(f"  {k}: AUROC={v['best_auroc']:.4f}  AUPRC={v['best_auprc']:.4f}  n_val={v['n_val']}")


if __name__ == "__main__":
    main()
