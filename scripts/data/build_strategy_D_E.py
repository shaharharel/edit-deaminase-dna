"""Build Strategy D (hard adversarial) and Strategy E (same-window).

Strategy E (same-protospacer-window proxy):
  For each positive at (chrom, pos), find OTHER Cs/As within ±20 bp on same strand
  that AREN'T in the positive set. These are positions in the same likely-guide window
  that could have been edited but weren't.
  - Same chrom, strand, ref base
  - Same gene (implicit since ±20 bp)
  - Different position
  - Trinuc may or may not match (we keep it loose — the protospacer geometry is the control)

Strategy D (hard adversarial):
  For each positive, train a quick XGBoost on (positives vs Strategy A negs) using
  the 81 handcraft features. Score ALL ±200 bp positions around each positive.
  Pick the highest-scoring NON-positive Cs as hard negatives.
  These are positions the simple model would mistake for OT.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import pyfaidx

HG38 = "/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
ROOT = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed")
POS_PARQUET = ROOT / "dna_offtarget_sites_v2.parquet"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def build_strategy_E(pos_df, fa, window=20, seed=42):
    """Same-window proxy: ±20 bp around each positive, find other valid Cs/As."""
    rng = np.random.default_rng(seed)
    pos_keys = set(zip(pos_df["chrom"], pos_df["pos"].astype(int)))
    negs = []
    failures = 0
    print(f"[E] sampling same-window (±{window}bp) negs for {len(pos_df)} positives...")
    for idx, p in pos_df.iterrows():
        chrom = p["chrom"]; pos = int(p["pos"]); strand = p["strand"]
        ref = p["ref"]
        # Look for ref bases in ±window around pos
        try:
            window_seq = fa[chrom][pos - 1 - window:pos - 1 + window + 1].seq.upper()
        except Exception:
            failures += 1; continue
        # For strand=-, we need bases on - strand: which means complement of + strand at the position.
        # Easier: scan the + strand window for the strand-aware ref base.
        # Strand-aware: if strand=+, look for ref. If strand=-, look for COMP[ref].
        target_genomic = COMP[ref] if strand == "-" else ref
        # Find all positions where window_seq[i] == target_genomic
        candidates = []
        for i, b in enumerate(window_seq):
            if b == target_genomic:
                cand_pos = pos - window + i
                if cand_pos == pos: continue  # skip self
                if (chrom, cand_pos) in pos_keys: continue  # skip other positives
                candidates.append(cand_pos)
        if not candidates:
            failures += 1; continue
        # Sample one random candidate from the window
        cand_pos = int(rng.choice(candidates))
        negs.append({
            "chrom": chrom, "pos": cand_pos, "strand": strand,
            "ref": ref, "alt": p["alt"],
            "editor": p["editor"], "source_paper": "matched_neg_E",
            "edit_rate": 0.0, "vaf": 0.0, "n_reps": 0,
            "source_tier": 1,
        })
        if (idx + 1) % 5000 == 0:
            print(f"  pos {idx+1}: negs={len(negs)} failures={failures}")
    print(f"[E] DONE. negatives={len(negs)} failures={failures}")
    return pd.DataFrame(negs)


def build_strategy_D(pos_df, fa, hand_features_path, seed=42, ratio=1.0):
    """Hard adversarial: use XGBoost on handcraft features to find high-scoring
    non-positive Cs in the ±500 bp window around each positive."""
    from sklearn.linear_model import LogisticRegression
    rng = np.random.default_rng(seed)
    pos_keys = set(zip(pos_df["chrom"], pos_df["pos"].astype(int)))
    print(f"[D] training quick LR model on handcraft features...")
    h = pd.read_parquet(hand_features_path)
    hcols = [c for c in h.columns
              if c not in {"chrom", "pos", "strand", "region", "center_trinuc"}
              and pd.api.types.is_numeric_dtype(h[c])]
    # Mark each handcraft row as positive if it's in pos_keys (vectorized)
    h_keys = list(zip(h["chrom"].values, h["pos"].astype(int).values))
    h["_is_pos"] = [k in pos_keys for k in h_keys]
    X = h[hcols].fillna(0).values
    y = h["_is_pos"].astype(int).values
    print(f"[D] LR train: {y.sum()} pos / {(1-y).sum()} neg")
    lr = LogisticRegression(max_iter=500, C=1.0, class_weight="balanced", n_jobs=-1)
    lr.fit(X, y)
    h["_score"] = lr.predict_proba(X)[:, 1]
    # For each positive, find high-scoring non-positives nearby (±10 kb)
    h_by_chrom = {c: g.sort_values("pos").reset_index(drop=True) for c, g in h.groupby("chrom")}
    negs = []
    failures = 0
    target_n = int(len(pos_df) * ratio)
    print(f"[D] sampling {target_n} hard adversarial negatives...")
    for idx, p in pos_df.iterrows():
        if len(negs) >= target_n: break
        chrom = p["chrom"]; pos = int(p["pos"])
        sub = h_by_chrom.get(chrom)
        if sub is None: failures += 1; continue
        nearby = sub[(sub["pos"] >= pos - 10000) & (sub["pos"] <= pos + 10000) & (~sub["_is_pos"])]
        if len(nearby) == 0: failures += 1; continue
        # Top-1 by LR score (most "positive-like" non-positive)
        best = nearby.nlargest(1, "_score").iloc[0]
        negs.append({
            "chrom": chrom, "pos": int(best["pos"]), "strand": best["strand"],
            "ref": p["ref"], "alt": p["alt"],
            "editor": p["editor"], "source_paper": "matched_neg_D",
            "edit_rate": 0.0, "vaf": 0.0, "n_reps": 0,
            "source_tier": 1,
            "lr_score": float(best["_score"]),
        })
        if (idx + 1) % 5000 == 0:
            print(f"  pos {idx+1}: negs={len(negs)} failures={failures}")
    print(f"[D] DONE. negatives={len(negs)} failures={failures}")
    return pd.DataFrame(negs)


def main():
    fa = pyfaidx.Fasta(HG38)
    pos = pd.read_parquet(POS_PARQUET)
    pos = pos[pos["chrom"].str.match(r"^chr(\d+|[XYM])$")].reset_index(drop=True)
    print(f"[main] {len(pos)} positives")

    # Strategy E
    neg_E = build_strategy_E(pos, fa)
    pos_out = pos.copy(); pos_out["is_positive"] = True
    neg_E["is_positive"] = False
    combined_E = pd.concat([pos_out, neg_E], ignore_index=True)
    combined_E.to_parquet(ROOT / "dna_offtarget_combined_E.parquet", index=False)
    print(f"[main] E: {len(combined_E)} total")

    # Strategy D (depends on handcraft features being downloadable from ai-chem)
    # We'll pull them via gcloud first
    hand_path = ROOT / "handcrafted_features_v3.parquet"
    if hand_path.exists():
        neg_D = build_strategy_D(pos, fa, hand_path)
        neg_D["is_positive"] = False
        combined_D = pd.concat([pos_out, neg_D], ignore_index=True)
        combined_D.to_parquet(ROOT / "dna_offtarget_combined_D.parquet", index=False)
        print(f"[main] D: {len(combined_D)} total")
    else:
        print(f"[main] handcraft features not local — skip D for now")


if __name__ == "__main__":
    main()
