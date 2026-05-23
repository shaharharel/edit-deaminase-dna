# DeaminaFormer-DNA — End-to-end guide-independent off-target prediction (2026-05-23)

## The validated pipeline
sources -> motif-QC'd canonical sites -> 1Mb bins + chromatin tracks -> accessibility model -> per-editor regional risk map.

## What we can predict (validated)
**The >=1 Mb guide-independent BE DNA off-target LANDSCAPE is predictable from chromatin accessibility:**
- LOCO Spearman **0.66** for both CBE (Yu+Lei) and ABE (Richter).
- Chromatin delta over opportunity+motif+mappability: **+0.09** [95% CI 0.07-0.11].
- **Cross-editor-class generalization: CBE-model->ABE 0.647, ABE-model->CBE 0.644** (~within-class).
- Controls all pass: circular-shift null ~0; leave-one-source-out 0.37-0.56; exceeds inter-assay ceiling 0.26;
  cell-type-robust 0.995 (HEK293/iPSC/HSPC maps near-identical).
- Driven by DNase/ATAC (open chromatin) + R-loop (ssDNA exposure) — mechanistically coherent.

## What we cannot (honest scope)
- Per-SITE / per-GENE which-base is stochastic (motif only; sub-1Mb non-reproducible).
- Absolute calibrated rates (no clean continuous DNA gold; needs bulk WGS).
- The deliverable is a >=1Mb RELATIVE-RANK screening prior (flags regions for deep sequencing; cannot clear a region).

## Deliverable
Per-editor regional off-target risk map (motif x accessibility), cell-type-robust, flagging cancer-driver
megabase regions (e.g. BE4: STK11, ERBB3, MEN1, TSC2 regions). data/processed/be4_risk_map.parquet.

## Bottleneck & next lever
DATA, not model (N=2,838 bins, label noise 0.26; deep nets/Enformer don't help, confirmed). Raw in-hand
bulk guide-independent sources are exhausted (Yu/Lei/Richter clean; others clonal/orthogonal/efficiency/mouse).
More accuracy requires DOWNLOADING new bulk genome-wide off-target WGS (Buchumenski studies; accessions in
their Supp S1) -> each adds a landscape, raises the 0.26 ceiling. + NB-GLM for calibrated uncertainty.
