# DeaminaFormer-DNA — Model Card / Honest-Claims Summary (v1, 2026-05-22)

Status: research prototype. NOT a release safety gate. Audit trail of what is and isn't validated.

## Objective
Predict **guide-independent** base-editor DNA off-target activity (intrinsic deaminase promiscuity on
exposed ssDNA) — the DNA analog of Levanon's RNA editing-index hotspots, but **predicted** (from
sequence + chromatin) rather than **measured**, on the permanent/heritable DNA substrate.

## Data
- Positives: ~70K published off-target sites; **guide-independent CBE** workhorses Doman + McGrath
  (rAPOBEC1), Yu (engineered), Chen (tdCBE); **guide-independent ABE** = Richter only.
- **Excluded:** CHANGE-seq (guide-dependent); **Lei** (coordinate bug — reads CCC not TpC; fix = −2 shift).
- Features: f = 81 handcrafted sequence features (trinuc/k-mer/GC/DNAshape); g = ENCODE accessibility
  (DNase, R-loop/MapR, ATAC, mappability). Gold = BE4_clone1 vs Parent_WGS treated/control index.

## What IS validated (defensible claims)
1. **Deaminase motif preference is reproducible & transferable.** rAPOBEC1 trinucleotide preference
   cosine **0.996** across cell types (HEK293 vs iPSC); cross-cell-type transfer AUROC **0.82** ≈ ceiling.
   Preference is **editor-specific** (engineered/dual deaminases differ → requires deaminase embedding).
2. **Guide-independent location reproduces only at ≥1 Mb.** Cross-assay Spearman: −0.25 (10 kb),
   −0.15 (100 kb), **+0.05 (1 Mb, crossover)**, +0.15 (5 Mb), +0.20 (10 Mb), chrom-block-bootstrap CIs.
   → Fine-scale (gene) location is NOT cross-assay reproducible; megabase domains are.

## What is NOT validated / honest limitations
- **Per-gene calibrated absolute burden is NOT cross-assay defensible** (gene-level Spearman negative).
  Per-gene output is only a **sequence-intrinsic prior** (accessibility-agnostic), explicitly flagged.
- **g (accessibility) is unidentifiable from one cell type** — gold is one editor (BE4), one clone,
  one cell type (HEK293). Cross-editor / cross-cell-type burden generalization is assumed, not shown.
- **ABE = single guide-independent source** → hypothesis-generating only; no transfer evidence.
  (Clinically inverted: ABE dominates in vivo therapy, yet has the weakest evidence here.)
- pos-vs-synthetic-negative AUROC is **confound-prone** (mappability/composition/source); replaced by
  gold-index correlation + cross-assay concordance as the validation backbone.

## Honest deliverable structure
1. Sequence-**intrinsic** per-gene susceptibility prior (reproducible, accessibility-agnostic).
2. **≥1 Mb** regional risk maps (the calibrated unit), cell-type-conditioned via g tracks.
3. Per-editor total-burden ranking (Buchumenski-style) — pending per-editor gold.

## QA-identified risks + mitigations (from ml-code-reviewer + scientific-analyst)
- Gold germline/clonal-SNP contamination → VAF filter + treated−control + **TpC-spectrum QC gate**.
- Negative leakage (random-genomic mappability/composition) → restrict to callable genome, trinuc+GC
  matching, permutation null, dinuc-shuffle (TODO for next model rev).
- Lei coordinate bug → drop now, re-extract with −2 shift later.
- Burden correlation → chromosome-block holdout + block-bootstrap CI + GC partial-correlation control.

## Validation backbone (in progress)
- P3 gold: predicted intrinsic burden (19,357 genes) vs BE4 WGS CDS C→T index — correlation pending
  (mpileup running; TpC QC gate before trusting the gold).
- Cross-substrate: Levanon RNA editing index (same deaminase) — pending RNA pipeline.

## Reproduce
`experiments/phase0_baseline/preference_transfer.py` (motif transfer);
`experiments/phase3_validate/resolution_concordance.py` (resolution curve);
`experiments/phase3_validate/{pred_gene_burden,correlate_gold}.py` (gold validation);
gold pipeline: `build_cds_bed.py` + `gold_mpileup_par.sh` + `aggregate_gold.py`.

---
## Update (2026-05-23): full autonomous run — the f/g verdict

After training the full model (incl. HyenaDNA on V100) and pursuing the g landscape:

**f (deaminase preference) = a trinucleotide motif. No ML needed.**
- trinuc-only (16 features) = 0.826 transfer ≈ full ±10bp/HyenaDNA models. Wider context & FM embeddings
  do NOT help (HyenaDNA 0.56-0.64 << 0.82). The guide-independent preference is TpC, a lookup table.
- Per-gene intrinsic susceptibility is ~uniform (TpC ubiquitous) → cancer drivers NOT enriched → per-gene
  flagging is near-content-free.

**g (accessibility landscape) = the genuine ML signal (modest, real).**
- Chromatin predicts the >=1Mb guide-indep off-target landscape beyond opportunity+ascertainment
  (LOCO Spearman 0.351 → 0.419), driven by DNase (open chromatin) + R-loop (ssDNA exposure).
  Generalizes across held-out chromosomes. NON-motif, learnable.
- Bounded by label reproducibility (~0.26 at 1Mb). Repli-seq (key ssDNA driver) collinear with DNase.

## Final verdict (honest)
- The defensible, ML-worthy product is a **regional (>=1Mb) chromatin-driven off-target susceptibility
  model**, per editor — NOT a per-site or per-gene predictor.
- The deaminase preference (f) is a motif: report it as a per-editor lookup + index, not ML.
- Per-editor total-burden ranking (Buchumenski-style) is the cleanest safety deliverable.

## Biggest levers (all acquisition/scope-gated)
1. Better LABELS: bulk treated+control deep WGS (the gold we lack) — raises the ~0.26 ceiling.
2. Replication timing + histone (modest, collinear with DNase).
3. Guide-DEPENDENT extension: predict per-therapy-guide off-targets (homology × window × context ×
   accessibility) — the ML-hard, clinically critical problem (CHANGE-seq-BE + libraries), set aside here.
