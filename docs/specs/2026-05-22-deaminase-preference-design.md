# DeaminaFormer-DNA — Deaminase-Preference Site-Rate Model (Design)

Date: 2026-05-22
Status: approved-to-build (foundational run in progress)
Author: Shahar + Claude (rethink session)

## 1. Problem reframe

Levanon's work (Buchumenski 2021; 2025 hotspot manuscript; Methods-in-Enzymology
protocol) **measures** guide-independent off-target editing **on RNA** via the editing
index (edits ÷ coverage over a region, treated − control). It is detection, not
prediction, and it is RNA.

This project goes two steps further: **(1) DNA** (permanent, heritable, FDA-relevant)
and **(2) prediction** instead of measurement. Goal: a regulatory-grade safety signal —
*which base editor deposits guide-independent DNA mutations, how much, and where*.

We model **guide-independent** off-targets only (intrinsic deaminase promiscuity on
exposed ssDNA). Guide-dependent off-targets (CHANGE-seq etc.) are a separate, guide-
sequence-dependent problem and are excluded from the burden model.

## 2. Evidence base (diagnostics run this session)

All on the 70,110-site pooled positive set + hg38 + refGene.

1. **Only 6.8% of pooled sites fall in CDS** (4,772). Called-site tables are sparse
   exactly where the clinic cares → CDS burden cannot be densely labeled from tables
   alone; needs the BAM index.
2. **Gene-level cross-source concordance fails** (guide-indep CBE: Lei/Yu/Doman
   Spearman −0.15 to −0.60). "Which gene is hot" does not reproduce across assays.
3. **Megabase-level concordance holds** (Lei↔Chen ρ=+0.27 p=2e-20; Lei↔Yu ρ=+0.25).
   Large-scale susceptibility (replication/R-loop domains) is reproducible.
4. **Deaminase motif preference is reproducible and transferable.** rAPOBEC1 BE4 in
   McGrath (iPSC) vs Doman (HEK293): C-centered trinucleotide preference cosine **0.996**,
   textbook TpC (TCT/TCA/TCC). ABE/Richter shows canonical TpA. → the **separable,
   transferable factor**.
5. **Preference is editor-specific** (Yu engineered YE1, Chen dual tdCBE differ) →
   requires the deaminase embedding + per-class heads.
6. **Lei is a QC red flag** — BE4max is rAPOBEC1 but reads CCC, not TC → residual
   strand/coordinate bug. Every source must be motif-QC'd vs its known deaminase.
7. **Data per deaminase:** rAPOBEC1 36,144 sites / 11 sources; TadA 32,710 / 2
   (1 guide-indep); rare variants <510 each.

### Interpretation
Site rate factorizes:
```
rate(site) ≈ f_deaminase(local sequence | editor)  ×  g_accessibility(R-loop/repliseq/ATAC | cell type)
```
- `f` (deaminase preference): reproducible, editor-specific, site-learnable. **Learn this.**
- `g` (accessibility): cell-type/condition specific. **Condition on target-cell tracks.**

Gene-level concordance failed because it measured the entangled product `f×g`. The
**motif-density component of gene burden is sequence-intrinsic and reproducible** even
across cell types; the absolute burden is `intrinsic × accessibility` and is validated
against ground truth (BE4 WGS index).

This also explains prior results: E=0.50 (no local determinism), AUROC≈0.65 (weak
site-vs-synthetic-negative), all consistent with a megabase-scale + motif-factor signal
that the old site-classification framing never isolated.

## 3. Modeling design

**Unit of prediction = site.** Best representation, ~10^7 candidate bases, captures the
transferable deaminase preference. Deliverable produced by **aggregation** to gene/exon.

- **Scope (denominator):** CDS + coding exons. ~9.7M C / ~8.8M A. Matches Levanon
  (cross-substrate comparability), clinically meaningful, tractable.
- **Two classes:** CBE (edited C, C→T/G→A) and ABE (edited A, A→G/T→C). Shared backbone,
  separate output heads + denominators.
- **Output:** calibrated per-base **edit rate** (not just presence). Where `edit_rate`/VAF
  exists (~34K sites) train as rate-regression; elsewhere presence with rate priors.
- **Aggregation:** expected gene/exon burden = Σ rate over eligible bases = predicted
  editing index (same quantity Levanon measures).

### Inputs
- **Local sequence (f):** start with handcrafted context (tri/penta-nucleotide one-hot,
  ±k window) — directly the deaminase motif; add HyenaDNA / Evo-2 embeddings later.
- **Accessibility (g):** ENCODE tracks already on ai-chem — R-loop (MapR), ATAC, DNase,
  histone, conservation, mappability. **Replication timing (repliseq) MISSING — download
  (key ssDNA-exposure determinant).**
- **Editor/deaminase embedding:** reuse shared decomposed schema (deaminase_family +
  mutations + cas + …) from `src/data/editor_features.py` → extrapolation to rare variants.

### No synthetic negatives
The genomic denominator (all eligible C/A) provides rate-0 background. Evaluation is by
**burden correlation / enrichment**, never synthetic-negative AUROC. This removes the
Strategy-D leak class entirely.

## 4. Labels & ground truth

- **Training (relative, many editors):** guide-independent published site tables, motif-QC'd
  per source. Drop guide-dependent (CHANGE-seq) and Doman-orthogonal-assay spatial signal
  (use Doman for *motif* only, not genome-wide spatial).
- **Gold calibration + held-out test (1 editor):** BE4_clone1 vs Parent_WGS deep WGS on
  ai-chem → real CDS C→T editing index (treated − control), via RNAEditingIndexer/mpileup.
  Recover the crashed Mutect2 call (missing `.stats`).
- **Cross-substrate validation:** Levanon RNA BAMs on ai-chem → RNA editing index per
  deaminase; test DNA↔RNA preference consistency for matched enzymes.

## 5. Evaluation

1. **Cross-source transfer (preference):** train deaminase-family A on source X, test on
   source Y same family (Doman→McGrath). Positive control. + cross-enzyme negative control
   (rAPOBEC1→TadA should NOT transfer).
2. **Burden correlation vs BE4 WGS index** (gold), held-out genes.
3. **Cross-substrate vs Levanon RNA** preference consistency.
4. **Per-editor burden ranking** (Buchumenski-style scalar) — robust deliverable.
5. **Calibration:** Brier, ECE < 0.05; reliability curve.
6. **Tier-A** (COSMIC/ClinGen) gene recall of the *intrinsic susceptibility* score.

## 6. Honest claims / limitations

- **Reproducible & defensible:** per-editor deaminase preference; gene-level *intrinsic
  (motif) susceptibility*; per-editor total burden ranking.
- **Context-dependent (validated, not assumed):** absolute gene burden = intrinsic ×
  cell-type accessibility; only validated where we have ground truth (BE4 WGS).
- **Weak/uncheckable:** guide-independent ABE transfer (single source); rare deaminases
  (embedding extrapolation, flagged low-confidence).
- All source labels motif-QC'd; Lei requires coordinate fix before inclusion.

## 7. Phased plan

- **P0 (now):** Foundational preference run — per-deaminase motif profiles + sequence-only
  site model with cross-source / cross-enzyme transfer AUROC. Decides: is preference
  learnable & transferable? (Expected: yes for rAPOBEC1.)
- **P1:** Source motif-QC + fix Lei; assemble CDS+exon candidate base index; download repliseq.
- **P2:** Factorized site-rate model (seq f + accessibility g + editor), C/A heads, calibrated.
- **P3:** Compute BE4 WGS CDS index (gold); validate aggregated burden.
- **P4:** Aggregate → gene/exon; Tier-A/B/C/D tiering; per-editor ranking; model card.
- **Backlog:** Evo-2 backbone; more treated/control WGS pairs; joint RNA+DNA.
