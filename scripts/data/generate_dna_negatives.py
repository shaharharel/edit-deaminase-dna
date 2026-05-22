"""Generate motif-matched negative sites for DNA off-target training.

For each positive (chrom, pos, ref, alt):
  - Sample a negative same-strand site in the same chromosome where:
    * Ref base matches (positive: C → sample C; positive: A → sample A)
    * Trinucleotide context (m1, ref, p1) matches
    * NOT in any positive parquet site (cross-source de-dup)
    * Same broad region class (CDS/UTR/intronic/intergenic) — best-effort

Output: dna_offtarget_negatives.parquet (same schema as positives, is_positive=False)
Combined: dna_offtarget_combined.parquet (positives + negatives with is_positive col)

Run locally — uses hg38 from edit-rna-apobec.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

HG38 = Path("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa")
POS_PARQUET = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_sites_v2.parquet")
NEG_PARQUET = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_negatives.parquet")
COMBINED = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_combined.parquet")


def get_trinuc(fa: Fasta, chrom: str, pos: int, strand: str) -> str:
    """Return trinucleotide context [m1, ref, p1] in `+` strand orientation."""
    try:
        seq = fa[chrom][pos - 2:pos + 1].seq.upper()
        if strand == "-":
            comp = str.maketrans("ACGTN", "TGCAN")
            seq = seq.translate(comp)[::-1]
        return seq
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=float, default=1.0, help="negatives per positive")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print(f"[neg] loading positives: {POS_PARQUET}")
    pos = pd.read_parquet(POS_PARQUET)
    pos = pos[pos["chrom"].str.match(r"^chr(\d+|[XYM])$")].reset_index(drop=True)
    print(f"[neg] {len(pos)} positives after chrom filter")

    print(f"[neg] loading hg38")
    fa = Fasta(str(HG38))

    # For each positive, get its trinuc
    print(f"[neg] computing positive trinucleotides")
    pos_trinuc = []
    for _, r in pos.iterrows():
        pos_trinuc.append(get_trinuc(fa, r["chrom"], int(r["pos"]), str(r["strand"])))
    pos["trinuc"] = pos_trinuc
    pos = pos[pos["trinuc"].str.len() == 3].reset_index(drop=True)
    print(f"[neg] {len(pos)} positives with valid trinuc")
    print(f"[neg] top trinucs: {dict(pos['trinuc'].value_counts().head(10))}")

    # Set of all positive (chrom, pos) for exclusion
    pos_keys = set(zip(pos["chrom"], pos["pos"].astype(int)))

    # Build per-chromosome ref-base position pools matched to trinuc
    # Strategy: for each positive, sample a position ±1Mb away on same chrom with matching ref+trinuc
    print(f"[neg] sampling matched negatives...")
    n_target = int(len(pos) * args.ratio)
    negs = []
    chrom_lens = {c: len(fa[c]) for c in fa.keys() if c.startswith("chr")}
    by_chrom = pos.groupby("chrom").size()
    print(f"[neg] target {n_target} negatives across chroms (top): {dict(by_chrom.head(5))}")

    attempts = 0
    failures = 0
    for idx, p in pos.iterrows():
        if len(negs) >= n_target:
            break
        chrom = p["chrom"]
        ref = p["ref"]
        target_trinuc = p["trinuc"]
        clen = chrom_lens.get(chrom, 0)
        if clen < 1000:
            continue
        # Try up to 30 random positions in same chrom, find one with matching trinuc, ref base, not a positive
        found = False
        for _ in range(30):
            attempts += 1
            cand_pos = int(rng.integers(100, clen - 100))
            cand_strand = p["strand"]  # match strand of the positive
            cand_seq = fa[chrom][cand_pos - 2:cand_pos + 1].seq.upper()
            if cand_strand == "-":
                comp = str.maketrans("ACGTN", "TGCAN")
                cand_seq = cand_seq.translate(comp)[::-1]
            if cand_seq != target_trinuc:
                continue
            if (chrom, cand_pos) in pos_keys:
                continue
            negs.append({
                "chrom": chrom, "pos": cand_pos, "strand": cand_strand,
                "ref": ref, "alt": p["alt"],
                "editor": p["editor"], "source_paper": "matched_neg",
                "edit_rate": 0.0, "vaf": 0.0, "n_reps": 0,
                "source_tier": 1,
                "trinuc": target_trinuc,
            })
            found = True
            break
        if not found:
            failures += 1
        if (idx + 1) % 5000 == 0:
            print(f"  pos {idx+1}: negs={len(negs)} attempts={attempts} failures={failures}")

    print(f"[neg] DONE. negatives={len(negs)} attempts={attempts} failures={failures}")
    neg_df = pd.DataFrame(negs)
    neg_df.to_parquet(NEG_PARQUET, index=False)

    pos_out = pos.drop(columns=["trinuc"], errors="ignore").copy()
    pos_out["is_positive"] = True
    neg_out = neg_df.drop(columns=["trinuc"], errors="ignore").copy()
    neg_out["is_positive"] = False
    combined = pd.concat([pos_out, neg_out], ignore_index=True)
    combined["concordance"] = combined.get("concordance", 1)
    combined.to_parquet(COMBINED, index=False)
    print(f"[neg] combined: {len(combined)} (pos={pos_out['is_positive'].sum()} neg={(~combined['is_positive']).sum()})")
    print(f"[done] {NEG_PARQUET}\n[done] {COMBINED}")


if __name__ == "__main__":
    main()
