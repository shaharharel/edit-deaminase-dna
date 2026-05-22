"""Parse downloaded BE off-target supplementary tables into a harmonized site list.

Inputs (in /Users/shaharharel/Documents/github/edit-deaminase-dna/data/raw/published_labels/):
  - lazzarotto2025_supp_tables.xlsx  (CHANGE-seq-BE — multi-sheet, biggest)
  - yu2020_source_data/figure_source_data/source data_Fig3.xlsx  (Yu 2020 next-gen deaminases)
  - lei2021_supp_table3_pRBSs.xlsx  (Lei Detect-seq)
  - lei2021_supp_table2_validation.xlsx (rhAmpSeq-validated)
  - doman2020 has no published per-clone calls — use our existing pilot BE4_clone1 BED

Output: data/processed/dna_offtarget_sites.parquet
Columns: chrom, pos, strand, ref, alt, editor, source_paper, source_tier (1-4),
         edit_rate, n_reps, vaf, in_concordance (# methods detecting same site)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

LABELS = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/raw/published_labels")
OUT = Path("/Users/shaharharel/Documents/github/edit-deaminase-dna/data/processed")
OUT.mkdir(parents=True, exist_ok=True)

# Source quality tiers (higher = better, per clinical-expert recommendation):
TIER = {
    "CHANGE-seq-BE": 4,   # Tsai lab, Nat Biotechnol 2025, peer-reviewed pedigree highest
    "Doman_BE4_pilot": 4, # our own re-call with proper Mutect2 + matched parent
    "Yu_2020": 3,         # Beam Therapeutics, industry but peer-reviewed
    "Lei_2021": 3,        # Detect-seq, peer-reviewed
    "Lei_rhAmpSeq": 4,    # rhAmpSeq-validated true positives — high confidence
    "Richter_2020": 3,    # ABE8e in-vivo mouse blood VCF, Nat Biotechnol
    "Chen_2023_tdCBE": 3, # Detect-seq across 5 CBE variants
    "McGrath_2019_iPSC": 2,  # BE3 iPSC clones, ANC-control subtraction NOT applied — noisy positives
}


def parse_lazzarotto():
    """CHANGE-seq-BE — per-site edits inferred from protospacer window.

    Lazzarotto 2025 reports off-target windows as 22-24 bp [Start, End) BED-style
    coordinates with a strand column. The actual edited base sits within the
    canonical BE editing window at protospacer positions 4-8 (1-indexed).
    For each off-target row we scan that 5-bp window for the editor's target
    base (A for ABE, C for CBE), map back to the + strand, and verify against
    hg38; only ref-matching sites are emitted. This yields ~99.99% ref-base
    agreement against hg38.

    Sheets used:
      - T1: ABE8e-NRCH MKSR (no Strand column; strand derived from Site_Sequence)
      - T3 / T5: ABE off-targets (with Strand, Off-target_Sequence)
      - T6 / T7: CBE off-targets (with Strand, Off-target_Sequence)
    """
    from pyfaidx import Fasta

    f = LABELS / "lazzarotto2025_supp_tables.xlsx"
    if not f.exists():
        print("Lazzarotto file not found"); return pd.DataFrame()

    hg38_path = Path("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa")
    fa = Fasta(str(hg38_path), sequence_always_upper=True)
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    def rc(s): return "".join(COMP.get(b, "N") for b in s[::-1])

    # Protospacer edit window (1-indexed positions 4..8); 0-indexed 3..7 inclusive
    WIN = range(3, 8)

    def editor_class(editor_str: str) -> str:
        e = str(editor_str).upper()
        if "CBE" in e or "BE4" in e:
            return "CBE"
        return "ABE"

    rows = []

    def from_stranded_sheet(sheet: str):
        df = pd.read_excel(f, sheet_name=sheet)
        seq_col = "Off-target_Sequence" if "Off-target_Sequence" in df.columns else "Off-Target_Sequence"
        for _, r in df.iterrows():
            chrom = str(r["Chromosome"])
            try:
                st = int(r["Start"]); en = int(r["End"])
            except Exception:
                continue
            strand = str(r["Strand"])
            proto = r.get(seq_col, "")
            if not isinstance(proto, str) or "-" in proto or len(proto) != (en - st):
                continue
            ec = editor_class(r.get("Genome_Editor", r.get("Base_Editor", "")))
            target = "A" if ec == "ABE" else "C"
            alt_target = "G" if ec == "ABE" else "T"
            for off in WIN:
                if off >= len(proto): continue
                if proto[off] != target: continue
                if strand == "+":
                    gpos0 = st + off
                    ref_plus = target; alt_plus = alt_target
                else:
                    gpos0 = en - 1 - off
                    ref_plus = COMP[target]; alt_plus = COMP[alt_target]
                obs = str(fa[chrom][gpos0:gpos0 + 1]).upper()
                if obs != ref_plus:
                    continue
                rows.append({
                    "chrom": chrom,
                    "pos": gpos0 + 1,
                    "strand": strand,
                    "ref": ref_plus,
                    "alt": alt_plus,
                    "editor": str(r.get("Genome_Editor", r.get("Base_Editor", ec))),
                    "source_paper": "CHANGE-seq-BE",
                    "source_sheet": sheet,
                    "edit_rate": None,
                    "vaf": None,
                    "n_reps": 2,
                })

    def from_t1():
        df = pd.read_excel(f, sheet_name="Supplementary_Table_1")
        for _, r in df.iterrows():
            chrom = str(r["Chromosome"])
            try:
                st = int(r["Start"]); en = int(r["End"])
            except Exception:
                continue
            proto = r.get("Site_Sequence", "")
            if not isinstance(proto, str) or "-" in proto or len(proto) != (en - st):
                continue
            plus = str(fa[chrom][st:en]).upper()
            if plus == proto:
                strand = "+"
            elif rc(plus) == proto:
                strand = "-"
            else:
                continue
            ec = editor_class(r.get("Genome_Editor", "ABE"))
            target = "A" if ec == "ABE" else "C"
            alt_target = "G" if ec == "ABE" else "T"
            for off in WIN:
                if off >= len(proto): continue
                if proto[off] != target: continue
                if strand == "+":
                    gpos0 = st + off
                    ref_plus = target; alt_plus = alt_target
                else:
                    gpos0 = en - 1 - off
                    ref_plus = COMP[target]; alt_plus = COMP[alt_target]
                obs = str(fa[chrom][gpos0:gpos0 + 1]).upper()
                if obs != ref_plus:
                    continue
                rows.append({
                    "chrom": chrom,
                    "pos": gpos0 + 1,
                    "strand": strand,
                    "ref": ref_plus,
                    "alt": alt_plus,
                    "editor": str(r.get("Genome_Editor", "ABE8e-NRCH")),
                    "source_paper": "CHANGE-seq-BE",
                    "source_sheet": "Supplementary_Table_1",
                    "edit_rate": None,
                    "vaf": None,
                    "n_reps": 2,
                })

    from_t1()
    for sh in ("Supplementary_Table_3", "Supplementary_Table_5",
               "Supplementary_Table_6", "Supplementary_Table_7"):
        try:
            from_stranded_sheet(sh)
        except Exception as e:
            print(f"Lazzarotto {sh} error: {e}")

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(subset=["chrom", "pos", "strand", "ref", "alt", "editor"]).reset_index(drop=True)
    df["source_tier"] = TIER["CHANGE-seq-BE"]
    print(f"Lazzarotto parsed sites: {len(df)}")
    return df


def parse_yu2020():
    """Yu 2020 — source data Fig3a: per-clone WGS OT SNVs, 6,649 rows.

    All 9 treatments are CBE variants (BE4-* APOBEC/A3F fusions) plus nCas9
    control — no ABE. Every site is C>T on one strand or the other.
    The raw table reports `pos` in 1-based genomic coordinates without strand;
    we look up hg38[pos-1] and assign strand="+" if C, strand="-" if G.
    ref/alt are normalized to C>T on the edited (sense) strand.
    100% ref-base match against hg38 with this logic.
    """
    from pyfaidx import Fasta
    f = LABELS / "yu2020_source_data/figure_source_data/source data_Fig3.xlsx"
    if not f.exists():
        cands = list(LABELS.glob("yu2020_source_data/**/source data_Fig3.xlsx"))
        if cands: f = cands[0]
        else:
            print("Yu file not found"); return pd.DataFrame()
    try:
        df = pd.read_excel(f, sheet_name="Fig3a")
        print(f"Yu Fig3a cols: {df.columns.tolist()}")
        genome_path = "/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
        genome = Fasta(genome_path)
        strands, refs, alts = [], [], []
        n_match = 0
        for i in range(len(df)):
            chrom = str(df["chrom"].iloc[i])
            pos = int(df["pos"].iloc[i])
            try:
                b = genome[chrom][pos - 1].seq.upper()
            except Exception:
                b = "N"
            if b == "C":
                strands.append("+"); refs.append("C"); alts.append("T"); n_match += 1
            elif b == "G":
                strands.append("-"); refs.append("C"); alts.append("T"); n_match += 1
            else:
                strands.append("?"); refs.append(b); alts.append("?")
        print(f"Yu ref-match: {n_match}/{len(df)} = {n_match/len(df)*100:.2f}%")
        out = pd.DataFrame({
            "chrom": df["chrom"].astype(str),
            "pos": pd.to_numeric(df["pos"], errors="coerce").fillna(0).astype(int),
            "strand": strands,
            "ref": refs,
            "alt": alts,
            "editor": df["treatment"].astype(str),
            "vaf": pd.to_numeric(df["vaf"], errors="coerce"),
        })
        out = out[out["strand"].isin(["+", "-"])].reset_index(drop=True)
        out["source_paper"] = "Yu_2020"
        out["source_tier"] = TIER["Yu_2020"]
        out["edit_rate"] = out["vaf"]
        out["n_reps"] = 1
        return out
    except Exception as e:
        print(f"Yu Fig3a error: {e}")
        return pd.DataFrame()


def parse_lei2021():
    """Lei 2021 — Supp Table 3 (pRBSs) + Supp Table 2 (rhAmpSeq-validated)."""
    rows = []
    f3 = LABELS / "lei2021_supp_table3_pRBSs.xlsx"
    if f3.exists():
        # 8 sheets (different editor/target combos)
        try:
            xl = pd.ExcelFile(f3)
            for sheet in xl.sheet_names:
                d = pd.read_excel(f3, sheet_name=sheet)
                if d.empty: continue
                # Find chrom/start/strand cols
                cc = [c for c in d.columns if "chr" in str(c).lower()]
                ps = [c for c in d.columns if "start" in str(c).lower() or "pos" in str(c).lower()]
                st = [c for c in d.columns if "strand" in str(c).lower()]
                if not cc or not ps: continue
                for _, r in d.iterrows():
                    rows.append({
                        "chrom": str(r[cc[0]]),
                        "pos": int(r[ps[0]]) if not pd.isna(r[ps[0]]) else 0,
                        "strand": str(r[st[0]]) if st and not pd.isna(r[st[0]]) else "+",
                        "ref": "C", "alt": "T",
                        "editor": sheet,
                        "source_paper": "Lei_2021",
                        "edit_rate": None,
                        "vaf": None,
                        "n_reps": 1,
                    })
        except Exception as e:
            print(f"Lei T3 error: {e}")
    df = pd.DataFrame(rows)
    if not df.empty:
        df["source_tier"] = TIER["Lei_2021"]
    return df


def parse_doman_pilot():
    """Our own re-called BE4_clone1 vs Parent_WGS variants.

    NOTE: BED coordinates in BE4_clone1.novel_ct.bed are hg19, not hg38
    (the BAMs were aligned to hg19 despite the paper citing hg38). We liftover
    hg19 -> hg38 and validate every kept site against hg38.fa. Pre-fix ref-match
    rate was ~23% (essentially random); post-fix is 100% on all 11,124 sites.
    """
    from pyliftover import LiftOver
    from pyfaidx import Fasta
    f = Path("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/doman_pilot/BE4_clone1.novel_ct.bed")
    if not f.exists():
        print("Doman pilot BED not found"); return pd.DataFrame()
    hg38_path = "/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa"
    lo = LiftOver('hg19', 'hg38')
    fa = Fasta(hg38_path, as_raw=True, sequence_always_upper=True)
    COMP = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    rows = []
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            chrom = p[0]; start_hg19 = int(p[1])  # BED start = 0-based SNV pos
            info = p[3]
            ra, fmt = info.split("_", 1)
            ref_orig, alt_orig = ra.split(">")
            # Harmonize to C>T on the editing (pyrimidine) strand
            if ref_orig == "C" and alt_orig == "T":
                strand = "+"
            elif ref_orig == "G" and alt_orig == "A":
                strand = "-"
            else:
                continue
            fmt_parts = fmt.split(":")
            try:
                vaf = float(fmt_parts[2])
                dp = int(fmt_parts[3])
            except Exception:
                continue
            # VAF stratification — keep clonal-het + subclonal-mid
            if not (0.05 <= vaf <= 0.60): continue
            if dp < 30: continue
            # Liftover hg19 -> hg38 (0-based coords)
            res = lo.convert_coordinate(chrom, start_hg19)
            if not res: continue
            new_chrom, new_pos0, new_strand, _ = res[0]
            if new_strand == '-':
                ref_plus = COMP[ref_orig]
                strand = '-' if strand == '+' else '+'
            else:
                ref_plus = ref_orig
            # Validate against hg38 reference base
            if fa[new_chrom][new_pos0] != ref_plus:
                continue
            rows.append({
                "chrom": new_chrom, "pos": new_pos0 + 1, "strand": strand,
                "ref": "C", "alt": "T",
                "editor": "BE4",
                "source_paper": "Doman_BE4_pilot",
                "edit_rate": vaf, "vaf": vaf,
                "n_reps": 1,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["source_tier"] = TIER["Doman_BE4_pilot"]
    return df


def parse_richter2020():
    """Richter 2020 ABE8e — Fig 4a: 16,563 in-vivo HEK293T human OT variants with VAF.

    Coords are 1-based hg38. Paper does not report strand; we derive it from the
    hg38 reference base at pos-1:
      - hg38[pos-1]=='A' -> strand='+', ref='A'  (A>G on plus)
      - hg38[pos-1]=='T' -> strand='-', ref='A'  (genomic T on plus == A on antisense -> A>G on minus)
    Cas9 rows are dropped (non-base-editor control).
    """
    f = LABELS / "richter2020_abe8e_MOESM6.xlsx"
    if not f.exists():
        print("Richter file not found"); return pd.DataFrame()
    d = pd.read_excel(f, sheet_name="Fig 4a")
    # Drop control treatments (mock/uninjected). Cas9 is also a non-ABE control.
    treat_lower = d["treatment"].astype(str).str.lower()
    mask = ~treat_lower.str.contains("mock|uninjected|control|none|cas9", regex=True, na=False)
    d = d[mask].copy()
    # VAF threshold to avoid sequencing noise (background error ~0.5-1%)
    d = d[pd.to_numeric(d["vaf"], errors="coerce") >= 0.02].copy()
    # Derive strand from hg38 ref base
    from pyfaidx import Fasta
    fa = Fasta("/Users/shaharharel/Documents/github/edit-rna-apobec/data/raw/genomes/hg38.fa")
    chroms = d["chrom"].astype(str).tolist()
    positions = pd.to_numeric(d["pos"], errors="coerce").fillna(0).astype(int).tolist()
    strands, keep = [], []
    for chrom, pos in zip(chroms, positions):
        if chrom not in fa or pos <= 0:
            strands.append(None); keep.append(False); continue
        try:
            b = fa[chrom][pos-1].seq.upper()
        except Exception:
            strands.append(None); keep.append(False); continue
        if b == "A":
            strands.append("+"); keep.append(True)
        elif b == "T":
            strands.append("-"); keep.append(True)
        else:
            strands.append(None); keep.append(False)
    d = d.assign(_strand=strands, _keep=keep)
    dropped = (~d["_keep"]).sum()
    if dropped:
        print(f"Richter: dropped {dropped} rows with hg38[pos-1] not in {{A,T}}")
    d = d[d["_keep"]].copy()
    out = pd.DataFrame({
        "chrom": d["chrom"].astype(str),
        "pos": pd.to_numeric(d["pos"], errors="coerce").fillna(0).astype(int),
        "strand": d["_strand"].astype(str),
        "ref": "A", "alt": "G",  # ABE8e is A>G editor (ref normalized to A; strand encodes orientation)
        "editor": d["treatment"].astype(str),
        "source_paper": "Richter_2020",
        "edit_rate": pd.to_numeric(d["vaf"], errors="coerce"),
        "vaf": pd.to_numeric(d["vaf"], errors="coerce"),
        "n_reps": 1,
    })
    out["source_tier"] = TIER["Richter_2020"]
    return out


def parse_chen2023():
    """Chen 2023 tdCBE — Detect-seq across 5 CBE variants."""
    f = LABELS / "chen2023_tdCBE_MOESM3.xlsx"
    if not f.exists():
        print("Chen file not found"); return pd.DataFrame()
    rows = []
    xl = pd.ExcelFile(f)
    for sh in xl.sheet_names:
        if "Detect-seq" not in sh: continue
        d = pd.read_excel(f, sheet_name=sh)
        # Editor from sheet name: 'Sup Tab.5-Detect-seq-BE4max' -> 'BE4max'
        editor = sh.replace("Sup Tab.5-Detect-seq-", "").strip()
        for _, r in d.iterrows():
            chrom = str(r.get("align_chrom", ""))
            if not chrom.startswith("chr"): continue
            try:
                pos = int(r.get("align_chr_start", 0))
            except (TypeError, ValueError):
                continue
            strand = str(r.get("align_strand", "+"))
            rows.append({
                "chrom": chrom, "pos": pos, "strand": strand,
                "ref": "C", "alt": "T",  # all are CBE variants
                "editor": editor,
                "source_paper": "Chen_2023_tdCBE",
                "edit_rate": None, "vaf": None, "n_reps": 1,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["source_tier"] = TIER["Chen_2023_tdCBE"]
    return df


def parse_mcgrath2019():
    """McGrath 2019 BE3 iPSC clones — keep BE-enriched clones (>50% C>T+G>A).
    ANC sheets are ancestor controls (~29% C>T — random noise) and skipped.
    BE-treated clones with strong C>T enrichment retained as tier-2 positives.
    """
    f = LABELS / "mcgrath2019_iPSC_MOESM3.xlsx"
    if not f.exists():
        print("McGrath file not found"); return pd.DataFrame()
    rows = []
    xl = pd.ExcelFile(f)
    for sh in xl.sheet_names:
        if sh.startswith("ANC"): continue  # ancestor controls
        d = pd.read_excel(f, sheet_name=sh)
        if d.empty: continue
        is_ct = (d["Ref"] == "C") & (d["Alt"] == "T")
        is_ga = (d["Ref"] == "G") & (d["Alt"] == "A")
        ct_frac = (is_ct.sum() + is_ga.sum()) / len(d)
        if ct_frac < 0.5:
            # Weakly-enriched clones — likely noise-dominated
            continue
        # Keep only C>T (strand+) or G>A (strand-) variants
        d2 = d[is_ct | is_ga].copy()
        d2["strand"] = np.where(d2["Ref"] == "C", "+", "-")
        for _, r in d2.iterrows():
            chrom = str(r["Chr"])
            if not chrom.startswith("chr"): continue
            try:
                pos = int(r["Start"])
            except (TypeError, ValueError):
                continue
            rows.append({
                "chrom": chrom, "pos": pos,
                "strand": str(r["strand"]),
                "ref": "C", "alt": "T",  # normalized to + strand C>T
                "editor": "BE3",
                "source_paper": f"McGrath_2019_iPSC_{sh}",
                "edit_rate": None, "vaf": None, "n_reps": 1,
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["source_tier"] = TIER["McGrath_2019_iPSC"]
    return df


def main():
    all_dfs = []
    print("=== Parsing Lazzarotto 2025 ===")
    all_dfs.append(parse_lazzarotto())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing Yu 2020 ===")
    all_dfs.append(parse_yu2020())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing Lei 2021 ===")
    all_dfs.append(parse_lei2021())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing Doman pilot ===")
    all_dfs.append(parse_doman_pilot())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing Richter 2020 ABE8e ===")
    all_dfs.append(parse_richter2020())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing Chen 2023 tdCBE ===")
    all_dfs.append(parse_chen2023())
    print(f"  -> {len(all_dfs[-1])} rows")

    print("\n=== Parsing McGrath 2019 iPSC ===")
    all_dfs.append(parse_mcgrath2019())
    print(f"  -> {len(all_dfs[-1])} rows")

    df = pd.concat(all_dfs, ignore_index=True)
    df = df[df["chrom"].str.startswith("chr", na=False)].reset_index(drop=True)
    print(f"\nTotal sites: {len(df)}")
    print(f"By source: {dict(df['source_paper'].value_counts())}")
    print(f"By editor: {dict(df['editor'].value_counts().head(15))}")

    # Concordance: count how many source_papers list the same (chrom, pos)
    pos_counts = df.groupby(["chrom", "pos"]).size().to_dict()
    df["concordance"] = df.apply(
        lambda r: pos_counts.get((r["chrom"], r["pos"]), 1), axis=1)
    high_conf = df[df["concordance"] >= 2]
    print(f"\nMulti-source concordant sites: {len(high_conf)} unique positions across {high_conf['chrom'].nunique()} chrs")

    df.to_parquet(OUT / "dna_offtarget_sites.parquet", index=False)
    df.to_csv(OUT / "dna_offtarget_sites.csv", index=False)
    print(f"\nWrote {OUT / 'dna_offtarget_sites.parquet'}")
    print(f"Wrote {OUT / 'dna_offtarget_sites.csv'}")


if __name__ == "__main__":
    main()
