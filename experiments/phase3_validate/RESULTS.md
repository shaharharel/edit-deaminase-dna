# P2/P3 — Validation results (2026-05-22)

## P2 v1 (DIAGNOSTIC ONLY — synthetic-negative confound, per QA)
f×g MLP, Strategy-B negatives, LOSO(CBE)/CV(ABE):
- CBE LOSO mean AUROC ≈ 0.59 (f) / 0.59 (f+g). Doman/McGrath ~0.54 held out — because Strategy-B
  trinuc-matching DELETES the motif signal (vs P0 random-neg 0.82). Yu=0.94, ABE-CV=0.96 are
  source-confound artifacts (global negatives ≠ per-source positive distributions).
- **Conclusion: pos-vs-synthetic-negative AUROC is confound-prone. Validation pivoted to gold-index
  correlation + cross-assay concordance.** (QA: ml-code-reviewer C1/C3, scientific-analyst P0-2.)

## P3a — Resolution-vs-concordance (THE honest-claim result)
Cross-assay reproducibility of guide-independent CBE burden (Doman/McGrath/Yu/Chen; Lei excluded
for coord bug), mean pairwise Spearman with chromosome-block-bootstrap 95% CI:

| Resolution | Spearman | 95% CI |
|---|---|---|
| 10 kb | −0.252 | [−0.255, −0.249] |
| 100 kb | −0.148 | [−0.151, −0.144] |
| **1 Mb** | **+0.052** | **[+0.039, +0.068]** (crossover) |
| 5 Mb | +0.152 | [+0.126, +0.185] |
| 10 Mb | +0.195 | [+0.139, +0.246] |

**Finding:** guide-independent off-target location is *anti-correlated* across assays at gene/10–100 kb
scale and only becomes reproducible at **≥1 Mb** (chromatin/replication domains). Combined with the
deaminase-motif reproducibility (cosine 0.996), the picture is:
- **local deaminase chemistry (motif): reproducible & transferable**
- **fine-scale location (≤100 kb): NOT reproducible across assays**
- **megabase regional susceptibility: reproducible**

**Honest deliverables:** (1) sequence-intrinsic per-gene susceptibility *prior* (accessibility-agnostic);
(2) ≥1 Mb calibrated regional risk maps; (3) per-editor total-burden ranking. Per-gene *calibrated*
absolute burden is NOT cross-assay defensible.

## P3b — BE4 WGS gold index: REJECTED by TpC QC (germline/clonal contamination)
Computed BE4_clone1 vs Parent_WGS over 186k CDS intervals (16.1M C sites, global CDS C→T 0.10%).
**TpC QC gate FAILED: 24.9% TpC** (APOBEC needs >40%; Doman/McGrath called sites = 81–90%).
Diagnosis (`clean_gold.py`): of 61,706 de-novo C→T variants (clean in parent),
- **98.6% are clonal-fixed (median VAF = 1.0, homozygous)** → germline/lineage differences, not edits
- TpC = 25% in BOTH clonal-fixed and low-level fractions → **no APOBEC enrichment at any VAF**
- only ~206 clean (TpC + sub-clonal) edits survive.
**Conclusion:** a single treated-clone vs single-parent WGS pair is germline-dominated and CANNOT serve
as a continuous editing-index gold. (Doman's *published* sites are clean only because Doman used
multi-clone/multi-control filtering.) A valid DNA editing-INDEX gold needs BULK deep WGS (not available).

### Validation backbone (revised, honest)
1. ✅ Cross-cell-type motif transfer (Doman↔McGrath, AUROC 0.82) — held-out lab+cell-type, same enzyme.
2. ✅ Resolution-concordance (≥1 Mb reproducibility).
3. ✅ Cross-source motif reproducibility (cosine 0.996).
4. ❌ BE4 WGS continuous gold — not usable (germline-dominated).
5. 🟡 Cross-substrate Levanon RNA index (same deaminase) — remaining external check.
The predicted intrinsic burden (19,357 genes, +GC control) stands as a sequence-intrinsic prior;
it has no clean DNA continuous ground truth to calibrate against with current data.
