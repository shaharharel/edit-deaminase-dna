# Step G — Bulk-WGS gold acquisition plan (2026-05-23)

Goal: a CLEAN continuous DNA editing-index gold (the g/burden + absolute-rate validation we lack).
Hard constraint learned: **clonal WGS is germline-dominated** (BE4_clone1 = 98.6% VAF≈1.0, 25% TpC).
A valid continuous index needs **bulk (non-clonal) deep WGS, treated + matched control**, OR
multi-clone/multi-control filtering (à la Doman's published-site pipeline).

## Candidate datasets (priority order)

1. **Buchumenski 2021 DNA-seq studies (D, K, I, J)** — PROVEN to yield a detectable DNA editing index
   (the only published DNA editing-index successes). Accessions in their Supplemental Table S1 (not in
   main text). Tool: `github.com/a2iEditing/BEIndexer`. ACTION: pull Supp S1 → SRA/GEO accessions →
   download → run BEIndexer (DNA mode) on our compute node. Highest-confidence path.
2. **2025 "Sensitive, direct detection of non-coding off-target base-editor unwinding/editing in
   primary cells"** (biorxiv 2025.09.25.678665) — modern, primary cells, direct (non-clonal) detection.
   Potentially the best contemporary bulk off-target dataset. ACTION: read data-availability, assess
   bulk vs enrichment, get accession.
3. **Yu 2020** (Nat Commun s41467-020-15887-5) — HEK293 CBE-variant WGS treated+control. We already have
   the called sites; check if RAW bulk WGS (not just clones) is deposited → could run the index. Likely
   clonal (same caveat) — verify.
4. **CHANGE-seq-BE 2025** (Nat Biotech s41587-025-02948-7) — in vitro, guide-DEPENDENT. NOT suitable for
   guide-independent burden; useful only as a guide-dependent reference. Deprioritize.

## What "valid gold" requires (acceptance criteria)
- Bulk (population) treated cells, NOT a single clone.
- Matched control (untreated or dead-deaminase), same lineage.
- Depth sufficient for low-level VAF (Buchumenski detection limit ~1e-5 → very deep, or rely on the
  aggregate index over many sites rather than per-site calls).
- **TpC-spectrum QC gate must pass** (>40% TpC for CBE) before trusting it — same gate that rejected the clone.

## Fallback if no clean bulk dataset exists
- Generate bulk treated+control deep WGS in-house (new sequencing) — the true gold; out of current scope.
- OR validate only the f factor (BE-DICT-style on-target efficiency, multi-editor) — strong but not genome-wide.

## Status
Acquisition is research+download (gated on accessions). The modeling side (f-preference) is validated;
the g/burden side remains gold-limited regardless of model quality.
