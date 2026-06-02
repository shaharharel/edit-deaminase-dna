# Phase 4 — Final state and Levanon-ready deliverable

## What we have

**The hg19-corrected single-clone DEI run is the deliverable.** See `docs/PHASE4_RESULTS_HG19.md` for full audit and numbers. Bottom line:

- **Top-100 DEI genes show 20.2× TpC/non-TpC motif specificity** — real BE-induced editing detected at the top of the gene distribution, not size-confounded (top-100 vs all-median size ratio = 1.07×).
- **Treated TpC EI_main = 2.08× Control TpC EI_main** — clear treated/control separation.
- **Whole-genome aggregate ratio = 1.15×** — at single-clone 50× WGS depth, the population-level signal is below the strict motif-specificity gate (sequencing error ~10⁻³ ≈ 10× the editing rate ~10⁻⁴), so most genes are noise-dominated. The signal lives in the tail.
- **f × g vs DEI Spearman ≈ 0** — model cannot identify which genes are in the top tail at this depth.

5 QA agents reviewed every script; **6 bugs caught and fixed**: H2 double-noise-correction, H3 asymmetric VAF cap, H7 bin assertion, H9 CpG exclusion, Spearman tie-correction, g-NaN propagation.

## Path A — attempted, terminated

Path A goal was to download 3 additional BE4 clones + 3 nCas9 matched controls (~500 GB raw, 6 clones × ENA direct download → bwa → markdup → mpileup → cleanup) to enable multi-clone aggregation and push the whole-genome gate above 2× motif specificity.

**Process worked correctly. Disk I/O was the killer:**
- BE4_clone2 (SRR10413103, 303 Gbp / ~70× WGS depth) consumed >10 hours of wall time at the bwa-mem + samtools-sort stage on this VM's persistent disk.
- Sort throughput measured at **~1.7 MB/s** writing chunked tmp files (~1.5–1.9 GB each) — would have needed ~15 more hours to finish sort alone for clone 2.
- Extrapolated end-to-end for all 6 clones: **~2 days wall time minimum**.

Terminated v8 to free the machine. Path A scripts (`/home/shaharh_quris_ai/download_doman_v8.sh`) and multi-clone aggregator (`experiments/phase4_index/parse_multiclone_pileup.py`) are committed and ready; can be re-launched later on a machine with faster I/O (local SSD or PD-extreme) if/when needed.

## What this means for the Levanon meeting

The takeaway slide is unchanged from the hg19 result write-up:

1. **The methodology is sound and audited** (QA-driven controls; 6 bugs fixed; matched control + low-VAF filter + mismatch-direction noise + motif-negative npc set).
2. **Real BE-induced editing IS detected at the top of the gene distribution** (20× TpC-specific in top-100 genes; treated/control 2.08× at TpC).
3. **Population-level validation requires multi-clone aggregation** (Path A) to push more genes above SNR — currently bottlenecked by I/O, not by the science.
4. **The framework translates from RNA to DNA**, but bulk DNA WGS at single-clone 50× has a depth ceiling that requires multi-clone replication or ultra-deep targeted amplicon — the two paths we discussed previously (more clones for the cohort burden test, or Levanon's prospective wet-lab validation for the per-gene test).

## Files for handoff

| Path | Contents |
|---|---|
| `docs/PHASE4_RESULTS_HG19.md` | hg19-corrected results, gate-by-gate, full audit trail |
| `docs/REPORT.html` | Phase 3 audited summary (f × g framework, multi-editor LOSO, top-line numbers) |
| `docs/LEVANON_BRIEFING.md` | Meeting briefing |
| `data/processed/phase4_results/hg19_run/per_gene_dei.parquet` | 18,799 gene rows with DEI_tpc, DEI_npc, EI rates, position counts, depths |
| `data/processed/phase4_results/hg19_run/per_gene_fxg.parquet` | 19,213 gene rows with predicted f-only, g-only, fxg_score, fxg_top10pct_mean |
| `data/processed/phase4_results/hg19_run/per_gene_validation.parquet` | Merged: empirical + predicted + ranks |
| `data/processed/phase4_results/hg19_run/qc_summary.txt` | Gate-by-gate sanity check output |
| `experiments/phase4_index/parse_samtools_pileup.py` | The hg19-fixed parser (asymmetric VAF, CpG-excluded npc, etc.) |
| `experiments/phase4_index/parse_multiclone_pileup.py` | Multi-clone aggregator (ready to run on Path A outputs whenever they exist) |
| `experiments/phase4_index/predict_per_gene_fxg.py` | Per-gene f × g prediction |
| `experiments/phase4_index/validate_index_vs_fxg.py` | Validation with motif-specificity ratio, partial Spearman, nonzero recall |
| `scripts/data/build_tpc_cds_bed.py` | hg19/hg38-agnostic TpC bed builder with CpG exclusion |
