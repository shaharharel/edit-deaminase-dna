# Phase 6 — editing-index enrichment (a–d) matrix across 4 datasets

Question: does a ~0.6–0.7 AUROC gradient-boosted seq+structure MODEL RECOVER the
guide-independent editing-index enrichment (Σed/Σcov over top-K% ranked, treated
vs control) — i.e. can a predictive model choose a burden-capturing off-target
panel — vs the ORACLE ceiling that ranks by the observed answer?

## Verdict: NULL everywhere. No model-recoverable, index-concentrating,
editor-specific signal in any dataset, by within-dataset OR cross-dataset ranking.

| Dataset | Editor | Assay | Oracle-specific (top1%) | Best MODEL AUROC | MODEL index/burden | Verdict |
|---|---|---|---|---|---|---|
| Doman   | BE4/YE1 CBE | WGS       | ~2.3× | 0.569 | 1.4% burden @top1% | null |
| Lei     | APOBEC CBE  | Detect-seq| 37×   | 0.50 (within-TpC) | ~1× | null |
| McGrath | AncBE4max CBE| WGS (isogenic ±induction) | 1.4× | 0.664 | germline-CpG artifact | null |
| Selict  | ABE8e       | enrichment| 5.5×  | 0.561 | 1.1× | null |
| (c) transfer | Doman↔Lei CBE | — | — | 0.48–0.50 (chance) | ≤ random baseline | null |
| (d) pooled   | = transfer (2 real-CBE) | — | — | 0.48–0.50 | ≤ random | null |

## Key adjudications
- **McGrath (`mcgrath_panel2.py`, `mcgrath_panel3.py`)**: the only stream with model
  AUROC in the 0.6–0.7 target range (0.664) is a **germline-CpG-deamination artifact,
  not editor**. Edited-C calls have ZERO TpC editor-motif enrichment (28% vs 29%
  background) but 3.6× CpG enrichment; the model's top-1% ranked panel is 68% CpG.
  Mechanism: the isogenic ±induction >3×-control filter at cov≥8 rejects balanced
  high-AF germline hets but concentrates the unbalanced low-AF background tail, which
  is mechanistically CpG-skewed. Fails the pre-registered TpC/nonTpC ≥1.5× bar (0.97×).
  Independently confirmed by scientific-analyst QA. (chr1-only pileup → spatial
  block CV, not held-out-chromosome; verdict rests on the compositional TpC-null.)
- **Cross-transfer (`cross_transfer.py`)**: harness-VALIDATED null. Within-Doman
  held-out-chromosome positive control recovers signal (AUROC 0.573), but Doman→Lei
  transfer is chance (0.498) and its index-enrichment tracks the built-in RANDOM
  baseline (1.76× vs random 1.80×) → no shared APOBEC editor sequence-preference
  transfers across datasets/assays. The DeaminaFormer shared-enzyme-transfer premise
  has no purchase here — off-target sites are diffuse/stochastic, not sequence-determined.

## Files
- `mcgrath_panel2.py` — McGrath oracle + spatial-block-held-out GB model panel.
- `mcgrath_panel3.py` — decisive CpG-vs-TpC decomposition of the McGrath 0.664.
- `cross_transfer.py` — (c) cross-transfer + (d) pooled, Doman↔Lei, with within-dataset
  positive control + random baseline.
- `ai_chem_startup_guarded.sh` — idempotency-guarded node startup (MCGRATH_WGS_DONE)
  to stop redundant reprocessing on spot-restart of the shared compute node.

Full logs / large intermediates (per-site parquets, pileups) live on GCS
`gs://ai-temp/apobec-genome-cache/{scripts,poc24_results}/`, not in git.

Deliverable stance: **measure per cell-context, don't predict** — a predictive
sequence panel for diffuse APOBEC/ABE guide-independent off-targets is not achievable
on this data; anchor on the A3A-hairpin enzyme-class dichotomy + lagging-strand covariate.
