# Phase 4 — DNA Editing Index pipeline: end-to-end result on Doman BE4 / Parent WGS

> Run date: 2026-06-01 (ai-chem CPU node). Honest negative result, fully audited.

## TL;DR

We built and ran an end-to-end Levanon-style DNA editing-index pipeline on the single available
BE4_clone1 + Parent_WGS pair (Doman 2020). The pipeline produces clean per-gene rates, the QA-mandated
controls all work, and **the result is an honest negative**: at single-clone WGS 50× depth, BE-induced
guide-independent editing is **not detectable above clonal somatic noise** in the per-gene aggregate.

- **Motif-specificity gate fails:** DEI(TpC) ≈ DEI(non-TpC) [mean ratio ≈ 0.78–0.93, gene-fraction = 0.41]
- **f × g vs empirical DEI Spearman ≈ 0.001** (and the motif-negative control ≈ 0.005 — same)
- **Recall@top-10% ≈ 9%** (random = 10%)

This does **not** invalidate the f × g model itself (which was already validated end-to-end on the
published guide-independent site catalogs in Phase 3). It says: **one BE4 clone vs one Parent at 50×
WGS is below the sensitivity threshold required to validate a per-gene burden prior in DNA**. The
additional Doman clones (~500 GB, downloadable) are the remediation.

## What we ran (5 stages)

| Stage | Script | Output |
|---|---|---|
| 1. CDS TpC inventory (CpG excluded) | `scripts/data/build_tpc_cds_bed.py` | 4.08 M TpC sites, 11.6 M non-TpC C's across 19,357 genes |
| 2. Parse pre-existing bcftools-AD TSVs + DEI aggregation | `experiments/phase4_index/parse_bcftools_index.py` | 18,162 genes with paired DEI |
| 3. Per-gene f × g prediction | `experiments/phase4_index/predict_per_gene_fxg.py` | 19,224 genes |
| 4. Position-count DEI (alternative metric) | `experiments/phase4_index/dei_position_count.py` | 18,309 genes |
| 5. Validation: Spearman + recall@top-K + motif-neg control | `experiments/phase4_index/validate_index_vs_fxg.py` | `validation_rate.parquet`, `validation_count.parquet` |

## Bugs caught and fixed (QA loop worked)

A QA agent reviewed every script statically before runs; a separate QA pass on results caught two more bugs during validation. **Six bugs fixed:**

| ID | File | Bug | Fix |
|---|---|---|---|
| H2 | `parse_bcftools_index.py` | Stacked treated-control on top of within-sample noise floor → mathematically double-subtracts | DEI = treated_main − control_main; noise floor reported separately as diagnostic |
| H3 | `parse_bcftools_index.py` | VAF cap applied per-sample → asymmetric germline filtering | Shared whitelist: drop position from both samples if VAF > cap in either |
| H7 | `predict_per_gene_fxg.py` | Bin lookup assumed integer index but no assertion → silent NaN propagation if bin = coordinate | Added assertion `b['bin'].max() < 1000` |
| H9 | `build_tpc_cds_bed.py` | CpG sites included in non-TpC control set → spontaneous deamination contaminates motif-negative | Excluded CpG (on both strands) from non-TpC control |
| – | `validate_index_vs_fxg.py` | Manual Spearman used no-ties formula → returned 1.0 for everything | Rewrote as Pearson-of-ranks (proper tie handling) |
| – | `predict_per_gene_fxg.py` | g composite used `.mean(axis=1)` not `nanmean` → NaN in any track → g_only all-NaN | nanmean + fallback for fully-NaN bins |

H8 (Levanon's true mismatch-direction noise floor would require strand-of-read; we use the transversion-average as a conservative proxy) was acknowledged as a methodological caveat. Did not block the result.

## The empirical DEI (post-fix)

| Metric | TpC | non-TpC | Ratio TpC/non-TpC |
|---|---|---|---|
| Treated `EI_main` (rate) | 1.18 × 10⁻⁴ | 1.18 × 10⁻⁴ | **1.00** ⚠️ |
| Control `EI_main` (rate) | 1.05 × 10⁻⁴ | 1.02 × 10⁻⁴ | 1.03 |
| **DEI (treated − control, rate)** | **1.27 × 10⁻⁵** | **1.63 × 10⁻⁵** | **0.78** ⚠️ |

| Sanity gate | Value | Need | Pass? |
|---|---|---|---|
| [a] DEI_tpc mean > 0 | ✓ | > 0 | ✅ |
| [b] Mean(DEI_tpc) / Mean(DEI_npc) | 0.78 | ≫ 1 | ❌ |
| [c] Fraction of genes with DEI_tpc > DEI_npc | 0.41 | > 0.70 | ❌ |
| [d] Treated EI_main > Control EI_main at TpC (gene-level) | 0.28 | > 0.50 | ❌ |

**Interpretation:** the post-VAF-filter rates *are* low and germline-realistic (~10⁻⁴ = 0.01%, in the right range for editing). But TpC and non-TpC C's show the **same** rate, in **both** treated and control. The treated-vs-control differential isn't motif-specific — it's a small clone-vs-parent difference that's the same at every C regardless of context.

**This is the noise floor of clonal somatic variation, not BE editing.**

### Why this happens at 50× WGS

- Per-position sensitivity at VAF 0.01 with 50× = 0.5 expected variant reads (mostly invisible)
- Per gene (≈ 280 TpC sites), expected variant reads from editing (rate ≈ 10⁻⁴ per C) ≈ 1.4
- Per gene, expected variant reads from sequencing error (rate ≈ 10⁻³ per C) ≈ 14
- **Noise floor ≈ 10× signal floor** at this depth. Per-gene aggregate cannot resolve.

The 206 high-confidence TpC edits in the pre-existing `be4_clean_edits.parquet` (VAF 0.05–0.31) confirm
*some* editing happened, but 206 events scattered across 19k genes is too sparse to compute per-gene rates.

## f × g vs empirical DEI

| Comparison | Spearman | p-value |
|---|---|---|
| DEI_tpc vs **f × g mean** | **0.001** | n.s. |
| DEI_tpc vs g_only | 0.002 | n.s. |
| DEI_tpc vs f_only | −0.011 | 0.13 |
| DEI_tpc vs `n_tpc` (gene size confound) | 0.042 | 2 × 10⁻⁸ |
| DEI_**npc** vs f × g mean (motif-neg control) | −0.007 | n.s. |
| Recall@top-10% of high-DEI genes by f × g | **8.9 %** | random = 10 % |
| Recall@top-1% of high-DEI genes by f × g | 0.6 % | random = 1 % |

**Zero predictive correlation, zero enrichment.** This is consistent with the QC finding above: the
empirical DEI does not contain BE editing signal, so there is no editing signal for the f × g model to
correlate with.

The only non-zero correlation is the **gene-size confound** (n_tpc, Spearman 0.04, p ≈ 2 × 10⁻⁸):
slightly larger genes have slightly more edit-candidate position counts, which makes their DEI slightly
non-zero. This is exactly the residual confound the QA reviewer flagged in advance.

## What this tells us (for Levanon)

**Honest framing of the experiment's outcome:**

1. **The methodology is sound.** The pipeline computes the right metric, the controls (VAF cap shared whitelist, CpG exclusion, treated−control subtraction, motif-negative control, Spearman with tie correction, NaN handling) all work, the QA loop caught every bug before it produced misleading numbers.

2. **The data is below sensitivity.** A single BE4 clone (≈ 50× WGS) compared to a single Parent (also ≈ 50× WGS) cannot give a clean per-gene editing-index signal because at this sequencing depth the BE-induced editing rate (≈ 10⁻⁴ per C) is buried under (a) sequencing error noise (≈ 10⁻³) and (b) clone-specific spontaneous somatic mutations indistinguishable from editing at the rate-aggregate level.

3. **The f × g model itself is not refuted** — Phase 3 already validated it end-to-end on the published guide-independent site catalogs (3–4× enrichment in held-out treatments and novel deaminases). What's refuted is the specific path of validating it via the per-gene DEI on this one clone.

4. **The remediation is more data, not a different model.** Two paths:
   - **Multi-clone Doman** (~500 GB downloadable from CLAUDE.md table: BE4 ×3 + YE1-BE4 ×3 + A3A-BE ×3, all sharing the Parent control). Per-position editing rate would average across multiple clones, suppressing clone-specific somatic; per-gene rates would integrate ≈ 9 × more reads at the same coverage.
   - **Targeted ultra-deep amplicon** (Levanon's specialty) at predicted high-f × g and predicted low-f × g genes, ≥ 1000× depth → would directly resolve sub-percent VAF and test whether predicted-high regions actually edit more.

5. **The Levanon RNA index works because** mRNA-seq has 1000× + depth on transcripts; editing rates in ADAR/APOBEC1 active regions are > 1%; many sites per transcript aggregate; matched paired controls remove systemic bias. **None of those four conditions hold in WGS DNA from one clone at 50×.** The framework translates conceptually; the per-base sensitivity required to make it run does not translate from RNA to bulk DNA.

## Files

| Path | Contents |
|---|---|
| `data/processed/phase4_results/per_gene_dei.parquet` | 18,162 genes × (rate-based DEI per gene, treated/control × TpC/nonTpC) |
| `data/processed/phase4_results/per_gene_dei_count.parquet` | 18,309 genes × (position-count DEI, VAF in 0.02–0.30) |
| `data/processed/phase4_results/per_gene_fxg.parquet` | 19,224 genes × (f_only, g_only, fxg_score_mean) |
| `data/processed/phase4_results/qc_summary.txt` | Position counts, VAF-filter survival, sanity-gate results |
| `data/processed/phase4_results/validation_rate.parquet` | Merged per-gene table for the rate-vs-fxg comparison |
| `data/processed/phase4_results/validation_count.parquet` | Same for position-count |

## Next concrete step (decision point)

**To make the DEI validation work, we need more BAMs.** Two options:

- **Option A: Download Doman additional clones** (~500 GB to ai-chem). The pipeline as-is then runs across BE4 ×3 + YE1-BE4 ×3 + A3A-BE ×3 — and we can cross-validate the f × g model against three different deaminases' empirical per-gene burdens. This is the strongest possible validation we can do without going to wet lab.

- **Option B: Wet-lab targeted ultra-deep sequencing with Levanon.** Pick ≈ 20 predicted high-f × g + 20 predicted low-f × g + 20 matched-control genes, deep-amplicon at ≥ 1000×, single editor (rAPOBEC1/BE4). Confound-proof prospective validation. This was already in the Levanon briefing (Ask #1).

The pipeline is now ready for either path — re-running it on additional clones requires only changing the BAM filenames in the bcftools mpileup step.
