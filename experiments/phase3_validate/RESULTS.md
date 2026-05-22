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

## P3b — Predicted intrinsic burden vs BE4 WGS gold index
- pred_gene_burden_v1: rAPOBEC1 preference model scored over CDS C's of 19,357 genes.
- Gold BE4 index: computing (BE4_clone1 − Parent_WGS, CDS C→T, TpC-QC gated). Correlation pending.
