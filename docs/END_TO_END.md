# DeaminaFormer-DNA — Guide-independent off-target prediction: honest end-to-end summary (2026-05-23)

> This supersedes earlier, over-claimed framings. Four QA passes (2 code + 2 scientific) found a
> gene-density confound running through the headline numbers; the claims below are the corrected,
> control-validated versions.

## Pipeline
sources -> motif-QC'd canonical sites -> 1 Mb bins + chromatin tracks -> regional model -> per-editor risk prior.

## What is actually true (control-validated)
1. **f (deaminase preference) = a trinucleotide motif (TpC/TpA).** No ML headroom (trinuc-only = full model;
   HyenaDNA loses). Per-editor, a lookup.
2. **Guide-independent off-targets concentrate in gene-dense, accessible megabase regions.** The >=1 Mb
   landscape is predictable at LOCO Spearman ~0.63 — but this is **dominated by gene density + sequence
   opportunity** (deterministic; gene-density alone r=0.61).
3. **Chromatin accessibility adds a MODEST, real increment beyond gene density: +0.04 LOCO**
   (delta over a baseline that includes GC + gene-density + motif + mappability). Down from a naive +0.09
   when gene density is NOT controlled. Autocorrelation-controlled (circular-shift null ~0).
4. **Cross-source generalization (leave-one-source-out): 0.37-0.56** — the honest cross-generalization number.
5. **Cell-type-robust** (DNase HEK293/iPSC/HSPC maps corr 0.995) — because 1 Mb accessibility is cell-type-
   invariant; the map carries NO target-cell-type specificity.

## What is NOT supported (corrected over-claims)
- "Landscape predictability 0.67 / chromatin +0.09" — the +0.09 was ~half gene-density confound; true
  accessibility increment is **+0.04**.
- "Cross-editor-class generalization 0.65 = learned coupling" — circular; the two landscapes co-localize in
  gene-dense regions (raw corr 0.55) and logA was a mirror of logC. Honest: one regional prior fits both.
- "Exceeds inter-assay reproducibility 0.26" — apples-to-oranges (in-target vs cross-assay). Dropped.
- "Flags cancer-driver regions" — mostly gene density (partial risk-vs-cancer | density = +0.10).
- Per-site / per-gene / absolute-rate prediction — not supported (motif stochastic; no continuous gold).

## Honest deliverable
A **>=1 Mb relative-rank screening prior** ranking megabase regions by (motif x gene-dense-accessible
localization). Useful to **triage regions for targeted deep sequencing**; **cannot clear** a region; not
calibrated; not target-cell-type-specific; cancer-content is largely gene-density-driven.

## Bottleneck & honest next steps
- Bottleneck = DATA (N=2,838 bins, label noise 0.26; deep nets/Enformer don't help). Raw in-hand bulk
  guide-independent sources exhausted (Yu/Lei/Richter).
- To strengthen: (1) more bulk genome-wide off-target WGS (EndoV-seq, Selict-seq, Buchumenski) -> raise ceiling;
  (2) gene-density-conditioned partial analysis as the standard; (3) NB-GLM calibration; (4) persist the bin
  builders into the repo (currently /tmp - reproducibility gap flagged by QA).
- Honest contribution: the f/g decomposition; guide-independent off-targets are a gene-dense-accessible-
  megabase phenomenon with a small specific accessibility increment; and the data-quality finding (only bulk
  genome-wide assays capture a predictable landscape).
