# CLAUDE.md - edit-deaminase-dna (BE DNA off-target safety gate)

## Environment

**Use the `quris` conda environment for all Python operations.**

## Project

**DeaminaFormer-DNA**: BE DNA off-target prediction tool aimed at being a regulatory-grade **safety gate** for base editor therapeutics.

Sister project: `edit-deaminase` (RNA off-targets, Levanon collaboration). This project covers the more ambitious clinical-safety direction.

## Why DNA matters more clinically

- RNA off-targets are transient (mRNA half-life hours-days). Risk: neoantigens, transient protein misfolding
- DNA off-targets are **permanent and heritable**. Risk: oncogenic transformation, haploinsufficiency, clonal expansion in ex vivo therapies (BEAM-101 HSPCs), regulatory concern for in vivo dosing (VERVE-101)
- FDA's January 2025 draft guidance for CRISPR products focuses on DNA off-target characterization. RNA off-targets are not yet a release-spec concern

## Architecture

**Separate from the RNA model** — the experts converged on Architecture B:
- Sequence encoder: **Evo (DNA foundation model)** instead of RNA-FM
- Structure features: **R-loop maps (DRIP-seq), replication timing (Repli-seq), chromatin accessibility (ATAC-seq)** instead of Vienna RNA folding
- **Shared with RNA model**: enzyme/BE-config embedding (the deaminase identity is substrate-agnostic)
- Substrate token: distinguishes DNA contexts (R-loop, replication fork, transient ssDNA hairpin)
- Multi-task heads: per-enzyme + per-substrate

Phase D (joint RNA+DNA) is **backlog only**. Build the two models separately and prove they work independently first.

## Training data (HUMAN only for V1; mouse deferred)

| Dataset | Editor | Status | Size |
|---|---|---|---|
| Lei 2021 Detect-seq | BE4max | Have | tiny BEDs |
| Doman 2020 BE4_clone1 + Parent_WGS | BE4 | Have (remote) | 268 GB BAMs |
| **Doman 2020 additional clones** | BE4, YE1-BE4, A3A-BE × 3 ea | Need to download | ~500 GB |
| **CHANGE-seq-BE 2025** (Lazzarotto/Tsai) | ABE8e, eA3A-BE3 genome-wide | Need to download | ~100 GB |
| Yu 2020 (Beam Therapeutics) | YE1-BE4-FNLS variants | Need | ~100 GB |
| Verve VERVE-101 (Rothgangl 2021) | ABEmax in NHP | Restricted — via Levanon | — |

## Strategy

| Phase | Goal | ETA |
|---|---|---|
| **A**: Re-evaluate v4_cds on Lei + Doman with new clinical lens (Tier A/B/C/D recall, calibrated) | 1 week |
| **B**: DeaminaFormer-DNA training on Lei + Doman + CHANGE-seq-BE + Yu | 2-3 weeks |
| **C**: Tier-stratified output spec, calibration, audit trail, model card | 1 week |
| **D (deferred)**: Joint RNA+DNA shared enzyme embedding ablation | backlog |

## Output spec (clinical safety gate)

Tier-stratified gene-category risk:
- **Tier A (hard block)**: COSMIC Tier 1 cancer (581) + ClinGen HI Level 3 (~340)
- **Tier B (flag and quantify)**: DepMap common essentials (~1,800) + tissue-specific
- **Tier C (annotate)**: All other coding + UTR + cCRE
- **Tier D (informational)**: Intergenic

Metrics: Tier A gene-stratified calibrated recall (NOT global AUROC). Brier score, ECE < 0.05.

## Validation

- APOB Q2153X recovery (canonical APOBEC1 substrate; sister project already passes this at top 1.6%)
- CHANGE-seq-BE held-out sites
- ClinVar pathogenic enrichment
- Per-patient germline VCF perturbation

## Shared code with edit-deaminase

- `src/data/editor_features.py` (same 33-editor catalog, locked EDITOR_ID_TO_IDX)
- `src/data/splitters.py` (leave_one_prjna_out, stratified_negative_sample)
- `src/utils/seeding.py` (set_seed)

These are duplicated, not symlinked, for repo independence. If they need updating, sync manually.

## Compute

- **CPU compute node** (large data volume) — primary DNA FASTQ + alignment
- **GPU node(s)** (T4 / V100 class) — model training when available
- **Local laptop** — model code development, smoke tests, small experiments
