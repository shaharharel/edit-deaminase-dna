# DeaminaFormer-DNA: Experiments Backlog & Plan

Updated: May 20, 2026.

Tracks ALL experiments we want to run, organized by phase. Status: ✅ done / 🟡 in progress / ⬜ TODO / ❌ dropped.

---

## Phase A — Baselines & Data (CLOSE BEFORE MOVING ON)

| ID | Experiment | Status | Outcome |
|---|---|---|---|
| A.1 | Harmonize all 8 published BE-OT label sources into single parquet (53K positives) | ✅ | `data/processed/dna_offtarget_sites.parquet` |
| A.2 | Strategy A negatives (trinuc + chrom matched) | ✅ | 21K negs, 0.40:1 |
| A.3 | Strategy B negatives (trinuc + region class) | ✅ | 19K negs, 0.36:1 |
| A.4 | Strategy C negatives (same-gene control) | ✅ | 36K negs, 0.68:1 |
| A.5 | HyenaDNA-small union embeddings (120K positions) | ✅ | `hyenadna_union.pt` |
| A.6 | Union ENCODE features (DNase, R-loop, ATAC, mappability) | ✅ | `encode_features_union.parquet` |
| A.7 | V1 sequence-only ablation table (A→A, B→B, C→C, B→C, A→C) | ✅ | C→C 0.66 / B→C 0.37 (chromatin confound shown) |
| A.8 | V2 sequence + ENCODE + DNase substrate | ✅ | C→C 0.663 (modest lift over V1) |
| A.9 | DNase-stratified calibration check on V2 | ✅ | PASSES — not chromatin-confounded |
| A.10 | Score Levanon hotspots with v4 endogenous model (Strategy A baseline) | ✅ | CBE 0.708 / ABE 0.804 gene-level |

---

## Phase B — Task-Specific Features (THE BIG GAP — current sequence model is generic)

| ID | Experiment | Status | Rationale |
|---|---|---|---|
| **B.1** | Hand-crafted motif features (24-d trinucleotide one-hot at edit position) | ⬜ | Canonical APOBEC/BE signal we used for RNA — missing here |
| **B.2** | Pentanucleotide context (5-mer one-hot, ~100 d) | ⬜ | Extended local sequence preference |
| **B.3** | GC content at ±10/50/100/500 bp windows | ⬜ | Sequence composition feature |
| **B.4** | k-mer composition (3/4/5-mer frequencies in ±200 bp) | ⬜ | Sequence environment fingerprint |
| **B.5** | Region class one-hot (CDS / 5UTR / 3UTR / intron / intergenic / ncRNA) | ⬜ | We classified for Strategy B — wire as feature |
| **B.6** | Distance to nearest exon boundary / splice site | ⬜ | Splice-edit confound; high signal |
| **B.7** | Distance to nearest TSS | ⬜ | Transcription start proximity |
| **B.8** | Conservation: phyloP + phastCons at edit site | ⬜ | Cross-species 24% enrichment already shown |
| **B.9** | DNAshape features (HelT, MGW, ProT, Roll) | ⬜ | DNA structural geometry from sequence |
| **B.10** | R-loop formation potential (RLFS sequence-based) | ⬜ | Complement to MapR experimental track |
| **B.11** | G-quadruplex propensity (QGRS or pqsfinder) | ⬜ | Local structure that affects BE access |
| **B.12** | Nucleosome occupancy prediction (NuPoP) | ⬜ | Sequence-based chromatin proxy |
| **B.13** | Replication timing (Repli-seq) — track downloaded but folder empty | ⬜ | Need to source — UCSC ENCODE |
| **B.14** | Histone marks (H3K27ac, H3K4me3, CTCF) | ⬜ | Need to download ENCODE tracks |
| **B.15** | Z-DNA / B-Z transition propensity | ⬜ | Sequence-derived structural |
| **B.16** | Gene expression context (GTEx TPM at the gene) | ⬜ | Highly-expressed genes more BE-vulnerable |
| **B.17** | chromHMM 15-state per site | ⬜ | Multi-track chromatin annotation |
| **B.18** | dbSNP common variant overlap (boolean exclude) | ⬜ | QA — exclude likely germline |
| **B.19** | RepeatMasker overlap (Alu, LINE, SINE) | ⬜ | Repeat context often confounds editing |
| **B.20** | mRNA stability features at the site (for downstream RNA model) | ⬜ | Half-life predictors |

---

## Phase C — BE-Specific Features (require guide RNA annotations — most papers don't share these)

| ID | Experiment | Status | Difficulty |
|---|---|---|---|
| C.1 | PAM context at protospacer position | ⬜ | Easy IF guide annotation exists per site |
| C.2 | Position within protospacer window (4-8 canonical edit window) | ⬜ | Need guide annotations |
| C.3 | Protospacer GC% | ⬜ | Need guide |
| C.4 | Guide-target Tm (melting temperature) | ⬜ | Need guide |
| C.5 | Guide self-complementarity / hairpin | ⬜ | Need guide |
| C.6 | Mismatch count between guide and target | ⬜ | Need guide |
| C.7 | Deaminase-active site distance prediction (BE3 ~10nt downstream of cut) | ⬜ | Sequence + guide |
| C.8 | UGI count (uracil DNA glycosylase inhibitor) → affects C-to-T efficiency | ⬜ | Per-paper editor identity |
| C.9 | Nickase vs full Cas9 (BE3 nicks; ABE7 doesn't) | ⬜ | Per-editor flag |
| C.10 | Linker length (BE3 has XTEN, ABE has shorter) | ⬜ | Per-editor |

---

## Phase D — Foundation Model Ablations (which backbone wins?)

| ID | Experiment | Status | Notes |
|---|---|---|---|
| D.1 | HyenaDNA-small 256d (current V2) | ✅ | C→C 0.663 |
| D.2 | HyenaDNA-medium 47M params / 160k context | ⬜ | T4 fp32 OK; bigger receptive field |
| D.3 | DNABERT-2 117M (after config registration fix) | ⬜ | BPE tokens, multi-species |
| D.4 | **Evo-2 1B** as DNA backbone | 🟡 | A100-b install nearly done; THE main comparison |
| D.5 | Evo-2 7B | ⬜ | Heavier, may not fit on T4 even fp16 |
| D.6 | Nucleotide Transformer v2 500M (with pinned transformers 4.39) | ⬜ | Earlier failed on `find_pruneable_heads_and_indices` import |
| D.7 | ESM-1b on translated peptide (CDS positions only) | ⬜ | Protein-context for coding-region OTs |
| D.8 | Frozen-backbone vs fine-tuned head | ⬜ | Ablation on capacity sharing |

---

## Phase E — Data Expansion (V3 datasets)

| ID | Experiment | Status | Yield estimate |
|---|---|---|---|
| E.1 | Liu 2022 GUIDE-tag (Joung lab, unbiased in-cell OT) | ⬜ | +10-30K sites |
| E.2 | Topham 2022 ONE-seq | ⬜ | +5-15K |
| E.3 | Tycko 2023 (Bhatt lab, thousands of guides × OT) | ⬜ | +10-30K |
| E.4 | Joung 2023 OutKnocker-seq | ⬜ | +5K |
| E.5 | Larger negative pools via pre-computed candidate index | ⬜ | bump B from 0.36:1 to 1:10 ratio |
| E.6 | Hard adversarial negatives (v4-model-scored high-prob non-positives) | ⬜ | "model can't trivially win" |
| E.7 | OutKnocker-seq + GUIDE-seq aggregator databases | ⬜ | Cross-source consistency |
| E.8 | Levanon BE-OT RNA per-site labels (REDItools) | 🟡 | RNA-specific; ~100K sites pending |

---

## Phase F — Model Architecture (after features land)

| ID | Experiment | Status | Hypothesis |
|---|---|---|---|
| F.1 | Multi-task head: OT probability + on-target rate (transfer from target-window data) | ⬜ | Auxiliary regularization |
| F.2 | Multi-task head: OT + REDIportal ADAR | ⬜ | Cross-substrate transfer |
| F.3 | Substrate-conditional editor embedding (DNA vs RNA branches) | ⬜ | Differentiate BE on DNA vs RNA off-target |
| F.4 | Cross-attention over chromatin context window (not just point) | ⬜ | Wider receptive field |
| F.5 | Distance-encoded position attention (BE prefers edit window pos 4-8) | ⬜ | Inductive bias |
| F.6 | Tier-stratified loss (upweight Tier-A cancer-driver positives) | ⬜ | Safety gate priority |
| F.7 | Calibrated probability output (Platt scaling / temperature) | ⬜ | Clinical interpretability |

---

## Phase G — Validation (after every major model rev)

| ID | Experiment | Status | Metric |
|---|---|---|---|
| G.1 | C→C AUROC + AUPRC + recall@1%FPR (current standard) | ✅ V2 done | sequence-only baseline |
| G.2 | DNase-stratified calibration | ✅ V2 done | PASSES |
| G.3 | Tier-stratified recall (Tier A genes only) | ⬜ | Clinical metric |
| G.4 | Real-genomic-distribution recall (score all hg38 Cs/As) | ⬜ | What the safety gate sees |
| G.5 | Cross-editor generalization (train BE3, test ABE) | ⬜ | Chemistry transfer |
| G.6 | Cross-cell-type generalization (train HEK293, test HSPC) | ⬜ | Substrate transfer |
| G.7 | Per-source ablation (drop-one-paper) | ⬜ | Robustness |
| G.8 | Calibration curve + Brier score + ECE | ⬜ | Clinical calibration |
| G.9 | Levanon hotspot recovery (gene-level rank correlation) | ⬜ | External validation |
| G.10 | ClinVar pathogenic enrichment at high-score sites | ⬜ | Clinical utility |
| G.11 | Per-PRJNA leave-one-out | ⬜ | Source generalization |

---

## Phase H — Tier-Stratified Clinical Output (after model is solid)

| ID | Experiment | Status |
|---|---|---|
| H.1 | COSMIC Tier 1 gene flag (581 genes) — Tier A | ⬜ |
| H.2 | ClinGen HI Level 3 (~340 genes) — Tier A | ⬜ |
| H.3 | DepMap essential common (~1800 genes) — Tier B | ⬜ |
| H.4 | Tier-stratified output spec & calibration | ⬜ |
| H.5 | Model card with audit trail | ⬜ |

---

## Phase Z — Joint RNA+DNA model (LAST PHASE — only after both close)

| ID | Experiment | Status | Note |
|---|---|---|---|
| Z.1 | Re-embed both RNA endo v3 (12K) + RNA BE-OT (Levanon) + DNA OT (53K+) with **Evo-2** | ⬜ | Evo-2 is DNA+RNA multi-modal — perfect fit |
| Z.2 | Joint training: shared editor embedding + substrate token discriminates DNA/RNA | ⬜ | Substrate token = endogenous_mRNA / BE_overexpr_mRNA_* / DNA_R_loop / DNA_replication_fork |
| Z.3 | Cross-modal validation: train on DNA, predict RNA OT (or vice versa) | ⬜ | Generalization test |
| Z.4 | Ablation: DNA-only vs joint — does joint training help DNA AUROC? | ⬜ | The whole point of unifying |

---

## Currently Running / Blocked

- 🟡 Evo-2 install on the A100 GPU node — model weights downloading
- 🟡 Levanon STAR alignment continuing on the compute node (REDItools deferred until per-PRJNA strand inferred)
- ⏸️ PRJNA923001 dropped (mouse + amplicon-seq, not transcriptome RNA-seq)
- ⏸️ Joint RNA+DNA (Phase Z) deferred until D.4 (Evo-2) lands + E.1-E.3 expansion done

---

## Priority queue for next session

1. **B.1-B.5** (motif + GC + region class) — quick wins, big expected lift
2. **D.4** (Evo-2 backbone) — gated on install
3. **B.8-B.11** (conservation + R-loop + G-quadruplex) — sequence-derived structural
4. **E.1-E.3** (more positive datasets) — data growth
5. **G.4** (real-genomic-distribution recall) — the clinically meaningful number
6. **G.3** (Tier-A recall) — safety gate metric
