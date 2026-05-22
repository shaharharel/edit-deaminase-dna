"""Strategies B + C for motif-matched DNA off-target negatives.

Strategy B (region-aware):
  For each positive, sample another position on the same chrom with:
  - Same trinucleotide context (m1, ref, p1)
  - Same ref base + strand
  - Same broad region class (CDS / 5UTR / 3UTR / intron / intergenic)
  - NOT in any positive set

Strategy C (same-gene):
  For each positive in gene G:
  - Find another C (or A) in the same gene's transcript that wasn't edited
  - Same ref base + strand
  - Trinucleotide matching encouraged but not required (small genes)
  Output: per-strategy parquet

Region class derived from RefGene at /Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt

Outputs:
  dna_offtarget_combined_B.parquet
  dna_offtarget_combined_C.parquet
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

HG38 = Path("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa")
REFGENE = Path("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/refGene.txt")
POS_PARQUET = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_sites_v2.parquet")
OUT_DIR = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed")


def load_refgene():
    """Parse UCSC refGene.txt to (chrom, exon_starts, exon_ends, cds_start, cds_end, strand, name)."""
    df = pd.read_csv(REFGENE, sep="\t", header=None,
                      names=["bin", "name", "chrom", "strand", "txStart", "txEnd",
                             "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                             "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"])
    df["exonStarts"] = df["exonStarts"].str.rstrip(",").str.split(",").apply(
        lambda x: [int(i) for i in x if i])
    df["exonEnds"] = df["exonEnds"].str.rstrip(",").str.split(",").apply(
        lambda x: [int(i) for i in x if i])
    return df


def classify_region(refgene_chrom: pd.DataFrame, pos: int) -> str:
    """Return 'cds', 'utr5', 'utr3', 'intron', or 'intergenic' for a position."""
    # Quick scan: find any transcript overlapping pos
    hits = refgene_chrom[(refgene_chrom["txStart"] <= pos) & (refgene_chrom["txEnd"] >= pos)]
    if hits.empty:
        return "intergenic"
    for _, t in hits.iterrows():
        # In exon?
        in_exon = False
        for s, e in zip(t["exonStarts"], t["exonEnds"]):
            if s <= pos <= e:
                in_exon = True; break
        if not in_exon:
            continue
        if t["cdsStart"] == t["cdsEnd"]:
            return "noncoding_exon"
        if t["cdsStart"] <= pos <= t["cdsEnd"]:
            return "cds"
        if t["strand"] == "+":
            return "utr5" if pos < t["cdsStart"] else "utr3"
        return "utr3" if pos < t["cdsStart"] else "utr5"
    return "intron"


def get_trinuc(fa: Fasta, chrom: str, pos: int, strand: str) -> str:
    try:
        seq = fa[chrom][pos - 2:pos + 1].seq.upper()
        if strand == "-":
            comp = str.maketrans("ACGTN", "TGCAN")
            seq = seq.translate(comp)[::-1]
        return seq
    except Exception:
        return ""


def strategy_B_region_aware(pos_df, refgene, fa, ratio=10.0, seed=42):
    """Sample negatives matched on trinuc + region class."""
    rng = np.random.default_rng(seed)
    # Classify all positives first
    print(f"[B] classifying {len(pos_df)} positive regions...")
    rg_by_chrom = {c: g for c, g in refgene.groupby("chrom")}
    pos_regions = []
    for _, r in pos_df.iterrows():
        rg = rg_by_chrom.get(r["chrom"])
        if rg is None:
            pos_regions.append("intergenic")
        else:
            pos_regions.append(classify_region(rg, int(r["pos"])))
    pos_df = pos_df.copy()
    pos_df["region"] = pos_regions
    print(f"[B] positive region distribution: {dict(pos_df['region'].value_counts())}")

    pos_keys = set(zip(pos_df["chrom"], pos_df["pos"].astype(int)))
    chrom_lens = {c: len(fa[c]) for c in fa.keys() if c.startswith("chr")}
    target_n = int(len(pos_df) * ratio)
    print(f"[B] sampling {target_n} negatives...")

    negs = []
    failures = 0
    for idx, p in pos_df.iterrows():
        if len(negs) >= target_n:
            break
        chrom = p["chrom"]; ref = p["ref"]; strand = p["strand"]
        trinuc = p["trinuc"]; pos_region = p["region"]
        clen = chrom_lens.get(chrom, 0)
        rg_chrom = rg_by_chrom.get(chrom)
        if clen < 1000 or rg_chrom is None:
            continue
        # Try up to 100 candidates
        found = False
        for _ in range(100):
            cand_pos = int(rng.integers(100, clen - 100))
            if (chrom, cand_pos) in pos_keys: continue
            cand_seq = fa[chrom][cand_pos - 2:cand_pos + 1].seq.upper()
            if strand == "-":
                comp = str.maketrans("ACGTN", "TGCAN")
                cand_seq = cand_seq.translate(comp)[::-1]
            if cand_seq != trinuc: continue
            cand_region = classify_region(rg_chrom, cand_pos)
            if cand_region != pos_region: continue
            negs.append({
                "chrom": chrom, "pos": cand_pos, "strand": strand,
                "ref": ref, "alt": p["alt"],
                "editor": p["editor"], "source_paper": "matched_neg_B",
                "edit_rate": 0.0, "vaf": 0.0, "n_reps": 0,
                "source_tier": 1,
                "trinuc": trinuc, "region": cand_region,
            })
            found = True
            break
        if not found: failures += 1
        if (idx + 1) % 5000 == 0:
            print(f"  pos {idx+1}: negs={len(negs)} failures={failures}")
    print(f"[B] DONE. negatives={len(negs)} failures={failures}")
    return pd.DataFrame(negs), pos_df


def _dist_to_exon_boundary(transcript, pos):
    """Min |pos - exon_start| over all exon starts/ends. Caps at 1e6."""
    md = 1e6
    for s in transcript["exonStarts"]:
        d = abs(pos - s)
        if d < md: md = d
    for e in transcript["exonEnds"]:
        d = abs(pos - e)
        if d < md: md = d
    return md


def _region_at(transcript, pos):
    """Region class within a transcript (cds/utr5/utr3/intron/ncrna)."""
    in_exon = any(s <= pos <= e for s, e in zip(transcript["exonStarts"], transcript["exonEnds"]))
    if not in_exon:
        return "intron"
    if transcript["cdsStart"] == transcript["cdsEnd"]:
        return "ncrna"
    if transcript["cdsStart"] <= pos <= transcript["cdsEnd"]:
        return "cds"
    if transcript["strand"] == "+":
        return "utr5" if pos < transcript["cdsStart"] else "utr3"
    return "utr3" if pos < transcript["cdsStart"] else "utr5"


def strategy_C_same_gene(pos_df, refgene, fa, ratio=1.0, seed=42):
    """Strategy C: same-gene negs. NOW match trinuc + region class + dist_to_exon_boundary
    (within ±100 bp window) — eliminates the dist_exon_boundary leak found at V5.
    """
    rng = np.random.default_rng(seed)
    rg_by_chrom = {c: g for c, g in refgene.groupby("chrom")}
    pos_keys = set(zip(pos_df["chrom"], pos_df["pos"].astype(int)))
    negs = []
    failures = 0
    DIST_TOL = 100  # bp — match dist_exon_boundary within this tolerance
    print(f"[C] sampling same-gene negatives (matched: trinuc + region + dist_exon ±{DIST_TOL}bp) for {len(pos_df)} positives...")
    for idx, p in pos_df.iterrows():
        chrom = p["chrom"]; ref = p["ref"]; strand = p["strand"]
        rg = rg_by_chrom.get(chrom)
        if rg is None:
            failures += 1
            continue
        hits = rg[(rg["txStart"] <= p["pos"]) & (rg["txEnd"] >= p["pos"])]
        if hits.empty:
            failures += 1
            continue
        t = hits.iloc[0]
        exon_positions = []
        for s, e in zip(t["exonStarts"], t["exonEnds"]):
            exon_positions.extend(range(s, e + 1))
        if len(exon_positions) < 100:
            failures += 1
            continue
        # Positive's anchor properties
        pos_region = _region_at(t, int(p["pos"]))
        pos_dist_exon = _dist_to_exon_boundary(t, int(p["pos"]))
        target_trinuc = p.get("trinuc")
        if not target_trinuc or len(target_trinuc) != 3:
            try:
                pos_seq = fa[chrom][int(p["pos"]) - 2:int(p["pos"]) + 1].seq.upper()
                if strand == "-":
                    pos_seq = pos_seq.translate(str.maketrans("ACGT","TGCA"))[::-1]
                target_trinuc = pos_seq
            except Exception:
                target_trinuc = None
        found = False
        for _ in range(400):  # bumped — 4-way match is harder
            cand_pos_1b = int(rng.choice(exon_positions)) + 1
            if (chrom, cand_pos_1b) in pos_keys: continue
            try:
                cand_trinuc = fa[chrom][cand_pos_1b - 2:cand_pos_1b + 1].seq.upper()
            except Exception:
                continue
            if strand == "-":
                cand_trinuc = cand_trinuc.translate(str.maketrans("ACGT","TGCA"))[::-1]
            if target_trinuc and cand_trinuc != target_trinuc:
                continue
            cand_center = cand_trinuc[1] if len(cand_trinuc) >= 2 else ""
            if cand_center != ref:
                continue
            # NEW: match region class
            cand_region = _region_at(t, cand_pos_1b)
            if cand_region != pos_region:
                continue
            # NEW: match dist_to_exon_boundary within tolerance
            cand_dist_exon = _dist_to_exon_boundary(t, cand_pos_1b)
            if abs(cand_dist_exon - pos_dist_exon) > DIST_TOL:
                continue
            negs.append({
                "chrom": chrom, "pos": cand_pos_1b, "strand": strand,
                "ref": ref, "alt": p["alt"],
                "editor": p["editor"], "source_paper": "matched_neg_C",
                "edit_rate": 0.0, "vaf": 0.0, "n_reps": 0,
                "source_tier": 1,
                "gene": t["name2"], "trinuc": target_trinuc,
                "region": cand_region, "dist_exon_boundary": float(cand_dist_exon),
            })
            found = True
            break
        if not found: failures += 1
        if (idx + 1) % 5000 == 0:
            print(f"  pos {idx+1}: negs={len(negs)} failures={failures}")
    print(f"[C] DONE. negatives={len(negs)} failures={failures}")
    return pd.DataFrame(negs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strategy", choices=["B", "C", "both"], default="both")
    args = ap.parse_args()
    pos = pd.read_parquet(POS_PARQUET)
    pos = pos[pos["chrom"].str.match(r"^chr(\d+|[XYM])$")].reset_index(drop=True)
    print(f"[main] {len(pos)} positives")
    fa = Fasta(str(HG38))
    print(f"[main] computing positive trinucs")
    pos["trinuc"] = [get_trinuc(fa, r["chrom"], int(r["pos"]), r["strand"]) for _, r in pos.iterrows()]
    pos = pos[pos["trinuc"].str.len() == 3].reset_index(drop=True)
    print(f"[main] {len(pos)} positives with valid trinuc")
    print(f"[main] loading refGene...")
    refgene = load_refgene()
    refgene = refgene[refgene["chrom"].str.match(r"^chr(\d+|[XYM])$")].reset_index(drop=True)
    print(f"[main] refGene: {len(refgene)} transcripts")

    if args.strategy in ("B", "both"):
        neg_B, pos_with_region = strategy_B_region_aware(pos, refgene, fa)
        combined_B = pd.concat([
            pos_with_region.drop(columns=["trinuc"], errors="ignore").assign(is_positive=True),
            neg_B.drop(columns=["trinuc"], errors="ignore").assign(is_positive=False),
        ], ignore_index=True)
        combined_B.to_parquet(OUT_DIR / "dna_offtarget_combined_B.parquet", index=False)
        print(f"[main] B: {len(combined_B)} total")

    if args.strategy in ("C", "both"):
        neg_C = strategy_C_same_gene(pos, refgene, fa)
        combined_C = pd.concat([
            pos.drop(columns=["trinuc"], errors="ignore").assign(is_positive=True),
            neg_C.assign(is_positive=False),
        ], ignore_index=True)
        combined_C.to_parquet(OUT_DIR / "dna_offtarget_combined_C.parquet", index=False)
        print(f"[main] C: {len(combined_C)} total")


if __name__ == "__main__":
    main()
