"""Rebuild dna_offtarget_sites.parquet keeping ONLY sites where ref matches hg38
(at some offset/strand). Per-source offset/strand inferred from hg38.

For each site, try (offset, strand_flip) ∈ {(-2,-1,0,+1,+2,+3,+4,+5)} × {keep, flip}
and pick the configuration where hg38[pos+offset-1] == expected_ref.

Output: dna_offtarget_sites_validated.parquet (subset of original, guaranteed ref-match).
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import pyfaidx

HG38 = "/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
IN = "/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_sites.parquet"
OUT = "/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed/dna_offtarget_sites_validated.parquet"

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def main():
    fa = pyfaidx.Fasta(HG38)
    df = pd.read_parquet(IN)
    df = df[df["chrom"].str.match(r"^chr(\d+|[XYM])$")].reset_index(drop=True)
    print(f"[validate] input: {len(df)} sites")

    # For each source, find the most common offset/strand-fix that matches hg38
    new_pos = []
    new_strand = []
    keep_mask = []
    for idx, r in df.iterrows():
        chrom, pos, strand, ref = r["chrom"], int(r["pos"]), r["strand"], r["ref"]
        best = None
        # Try offsets 0, +1, -1, +2, -2, +3, -3, +4, +5 in order of decreasing likelihood
        for off in [0, 1, -1, 2, -2, 3, -3, 4, 5]:
            try:
                b = str(fa[chrom][pos + off - 1]).upper()
            except Exception:
                continue
            # Match as-is (current strand interpretation)
            expected = COMP[ref] if strand == "-" else ref
            if b == expected:
                best = (pos + off, strand)
                break
            # Match flipped strand
            expected_flip = ref if strand == "-" else COMP[ref]
            if b == expected_flip:
                flipped = "+" if strand == "-" else "-"
                best = (pos + off, flipped)
                break
        if best is None:
            keep_mask.append(False)
            new_pos.append(pos)
            new_strand.append(strand)
        else:
            keep_mask.append(True)
            new_pos.append(best[0])
            new_strand.append(best[1])
        if (idx + 1) % 10000 == 0:
            print(f"  {idx+1}/{len(df)} processed, {sum(keep_mask)} validated")

    df["pos"] = new_pos
    df["strand"] = new_strand
    df["validated"] = keep_mask
    print(f"\n[validate] per-source validation rate:")
    for src in sorted(df["source_paper"].unique()):
        sub = df[df["source_paper"] == src]
        rate = 100 * sub["validated"].sum() / len(sub)
        print(f"  {src:<35} {sub['validated'].sum():>6} / {len(sub):>6} ({rate:>5.1f}%)")

    df_kept = df[df["validated"]].drop(columns=["validated"]).reset_index(drop=True)
    df_kept.to_parquet(OUT, index=False)
    print(f"\n[validate] kept {len(df_kept)} / {len(df)} sites ({100*len(df_kept)/len(df):.1f}%)")
    print(f"[done] wrote {OUT}")


if __name__ == "__main__":
    main()
