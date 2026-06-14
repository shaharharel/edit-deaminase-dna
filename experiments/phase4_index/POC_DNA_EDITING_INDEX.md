# POC: DNA editing index for BE4 off-target (Doman 2020, hg19)

**Date:** 2026-06-14 · **Status:** measurement POC complete (predictive model deferred to Phase 2)

## Goal
Translate the Levanon RNA editing-index principles to **DNA** off-target editing and test, on
data in hand, whether a base editor's genomic off-target burden is **measurable and reproducible** —
before investing in the full cross-deaminase cohort or a predictive model.

## Data
Doman 2020 cross-deaminase WGS, CDS C/G candidate sites (hg19), pooled samtools mpileup.

| role | samples |
|---|---|
| treated (BE4 = rAPOBEC1) | BE4_clone1, BE4_clone5, BE4_clone6, BE4_clone7 |
| control (untreated) | Parent_WGS |

All 4 treated clones independently QA-passed (per-clone TpC C→T 1.81–2.46e-4 vs Parent ~2.1e-4).

## Method (Levanon-faithful)
- **Pooled (read-summed) rate estimator**, treated minus matched control.
- **No per-site coverage floor.** A site is kept if pooled-treated coverage ≥1 **and** control ≥1 —
  low-coverage sites are *not* dropped; they accumulate across clones. (This is the key change vs the
  prior `parse_multiclone_pileup.py`, which took the intersection across all samples and silently
  discarded any site missing in even one clone.)
- Germline whitelist: drop a site iff Parent VAF > 0.05 (inherited SNP, not low coverage).
- Motif channels: **TpC** (APOBEC1 signature) vs **non-TpC** (negative control / per-gene noise proxy).
- Gene filter for per-gene stats: ≥20 TpC positions.

Scripts: `poc_dna_editing_index.py` (main aggregation) · `poc_splithalf.py` (reproducibility).
Outputs in `poc_results/` (`per_gene_dei.parquet`, `per_site_labels.parquet`, `burden_capture.csv`,
`qc_summary.txt`, `splithalf.txt`, `figures/`).

## Findings (honest, mixed → net positive)

### 1. Aggregate signal is real
Pooled per-gene treated TpC rate **6.2e-5** vs control **2.7e-5**; 65.6% of genes show treated > control.
Consistent with per-clone QA (1.15–1.26× at the bulk read level).

### 2. The signal is extraordinarily sparse (the deep-Levanon regime)
Even pooling 4 clones at **322× mean depth**, only **1.45%** of TpC sites carry a single edited read
(mean 0.02 edits/site). ~98.5% of sites are zeros. This is exactly the accumulation regime — no single
site is informative; burden emerges only by pooling.

### 3. Per-gene index is reproducible *only with pooling* — and it is
- **Single-clone** per-gene rate reproducibility: Spearman ≈ 0.13 (too sparse — not usable alone).
- **Split-half** (pool clone1+5 vs clone6+7): **Spearman 0.41** (DEI_tpc), **0.47** (motif-corrected
  TpC−nonTpC), both p≈0. Top-gene overlap 3.7× enriched over random; half-B DEI rises monotonically
  across half-A quintiles (−4.9e-5 → +6.9e-5).
- **The motif-corrected index reproduces *better* than the raw index** → the reproducible component is
  the APOBEC-specific (TpC) part, exactly as expected.

→ Pooling is **essential**, not optional. With it, the per-gene DNA editing index is a real,
moderately-reproducible quantity. More clones/depth will raise 0.41 further.

### 4. Off-target burden is concentrated (targeted-panel claim holds)
50% of excess TpC edit burden is captured by the top **9.6%** of genes; 80% by 29%; 90% by 41%.
A small high-yield gene panel captures most of the burden — the Levanon→DNA "reduced search space /
precision@K" translation.

### 5. Caveat against the naive top-K specificity number
The headline "top-100 = 9.6× TpC-vs-nonTpC" is **selection on the pooled data** (winner's curse): the
*exact* top genes are only modestly stable across clones (Jaccard 3.7× random). The trustworthy claims
are the **split-half rank correlation (0.41–0.47)** and the **burden concentration**, not the precise
identity of any individual top gene at this depth.

## Figures
- `figures/splithalf_reproducibility.png` — independent-half per-gene rate agreement
- `figures/motif_specificity_topK.png` — TpC vs non-TpC DEI by top-K
- `figures/burden_capture.png` — precision@K burden-capture curve

## Implications for the plan
1. **Aggregate + gene-level editing index = robust deliverable now.** The DNA editing index reproduces
   at gene resolution when pooled; the motif-corrected version is the cleanest readout.
2. **Site-level localization is depth-limited.** Per-site prediction needs the count-likelihood model
   (Beta-Binomial: low-coverage sites contribute without thresholding) + more depth/clones. Evaluate at
   gene/aggregate level, not per-site accuracy.
3. **More clones are worth it** — each one tightens the per-gene correlation (0.13 single → 0.41 pooled-2
   → higher with more). We are not accumulating clones without return.
4. **Phase-2 model** (site+gene hierarchical, gene = differentiable sum of site head, Beta-Binomial loss)
   is justified at **gene** resolution; the per-site/per-gene label tables produced here
   (`per_site_labels.parquet`, `per_gene_dei.parquet`) are its training targets.
5. **Cross-deaminase contrast** (BE4 vs YE1-BE4 vs nCas9) remains the headline differentiator and is the
   reason to continue acquisition — none of those clones have landed yet.

## Open follow-on (not in this POC)
- **PTC hotspots** (CGA→TGA, CAA→TAA, CAG→TAG): needs hg19 CDS reading-frame + reference context, which
  the C-only mpileups don't carry. Requires a CDS-frame annotation join — queued as the next step.
