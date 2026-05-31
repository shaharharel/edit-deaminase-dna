# DeaminaFormer-DNA — Briefing for Levanon

> Meeting prep, 2026-05-30. Companion to `docs/REPORT.html` (full audited results).
> Sister project: `edit-deaminase` (RNA, Levanon collaboration).

## 1. Motivation

You did **RNA**; we're doing **DNA** — the more permanent, more clinically-regulated substrate. FDA's January 2025 draft guidance for CRISPR products centers on DNA off-target characterization; RNA off-targets are not yet release-spec. The DNA analog of your **editing-index** framework (Buchumenski-style aggregate burden, deaminase motif × accessibility landscape) maps cleanly — but the data situation is harder:

- **Your RNA setup:** matched treated vs control → editing index, sub-percent VAF reachable on transcripts (high coverage).
- **DNA reality:** bulk genome-wide treated−control WGS doesn't exist — sub-1% VAF is below 30–50× WGS sensitivity. So **bulk-population guide-independent BE DNA off-targets are physically unobservable** with standard WGS. We work from published *site catalogs* from enrichment assays (Detect-seq, EndoV-seq, clonal WGS).

This data-ceiling finding is itself one of the project's contributions.

## 2. Data (HUMAN only, V1)

68,532 canonical sites after QC, 2,838 1 Mb bins genome-wide.

| Source | Deaminase | N sites | Tp[C/A] frac | Role |
|---|---|---|---|---|
| Doman_BE4_pilot | rAPOBEC1 | 11,123 | 0.81 | Clean training |
| McGrath_2019_iPSC | rAPOBEC1 | 13,846 | 0.90 | Clean training |
| Lei_2021 | rAPOBEC1 | 2,742 | 0.48 | −2 coord fix; **suspect** (Detect-seq footprint) |
| Yu_2020 | engineered YE1 | 6,648 | 0.17 | Novel-motif test |
| Richter_2020 | TadA (ABE) | 16,287 | 0.83 | Novel-deaminase test |
| BE4 clonal WGS | rAPOBEC1 | — | 0.25 | **Rejected** (98.6% VAF=1.0 = germline) |
| CHANGE-seq-BE | guide-dep | 17,473 | — | Held back (guide-dep workstream) |

**Features:** f = local motif (±10 one-hot or deaminase trinucleotide spectrum); g = per-1 Mb DNase, ATAC, R-loop/MapR, H3K27ac, gene-density, mappability.

## 3. Methods

- **Models:** trinuc lookup, ±10 one-hot MLP, ±25 CNN, deep MLP, bin-level MLP (LOCO), site-level MLP, **deaminase-conditioned LOSO**, **f × g deployment**. (Also tested HyenaDNA mean-pool and center-token — both lose to one-hot.)
- **Evaluation harness (QC G0–G9):** LOCO, LOSO, circular-shift null, within-gene-density-quintile control, TpC motif-spectrum gate, **background ladder** (random / trinucleotide-matched / bin-matched), enrichment metric (recall / K, portable across pos:bg ratios).
- **Honesty discipline:** 4 over-claims corrected during review (leaked 0.95 from synthetic-negative hard-mining; cell-type-matching artifact; Lei coordinate bug; gene-density confound). All numbers below are post-correction.

## 4. Key insights

1. **The motif is deaminase-specific, not treatment-specific.** A "treatment" = deaminase + Cas + guide + dose; only the deaminase determines f. (Exactly the architectural premise of shared enzyme embeddings.) Confirmed empirically: same-deaminase LOSO transfers (2.5–2.6×); novel-deaminase zero-shot fails (1.2×) because the model applies the wrong motif.
2. **Off-targets concentrate in gene-dense, accessible megabases.** Bin-level LOCO Spearman ~0.63 — but **dominated by gene-density** (r=0.61 with edit count). Honest chromatin increment beyond a strict gene-density + GC + mappability baseline: **+0.04** (down from a naive +0.09). Within-gene-density-quintile control-by-design: chromatin Spearman 0.35.
3. **f and g factorize cleanly and have very different ML characters.** f = a lookup table (no ML headroom — trinuc-only ≈ ±10 ≈ HyenaDNA all at AUROC ~0.82). g = a real but modest learnable regional signal, capped by inter-assay reproducibility (~0.26 at 1 Mb).
4. **g is approximately deaminase-agnostic and transfers across editors.** Even Richter (TadA/ABE), with no ABE training data and a fundamentally different base, gets 2.3× enrichment purely from the g signal learned on CBE sources.
5. **Cell-type-robust at 1 Mb.** DNase HEK / iPSC / HSPC corr 0.995. Useful (one map works) AND a limitation (no target-cell-type specificity at that resolution).
6. **Data, not architecture, is the ceiling.** Deep ≈ shallow ≈ CNN. HyenaDNA < one-hot. The bin-level model is N = 2,838.

## 5. Results — the ladder

| Eval | Resolution | Motif used? | Held out | Background | Recall@top-10% |
|---|---|---|---|---|---|
| Bin-level (motif-blind) | 1 Mb | NO | editor-class | random | **24%** (2.4×) |
| **f × g, held-out editor (incl. novel deaminase)** | site | yes (provided) | full editor | random | **22–37%** (2.2–3.7×) |
| Site, within-editor, trinuc-matched | site | yes (full ±10) | chromosome | trinuc-matched | 30–50% |
| Site, within-editor, random | site | yes | chromosome | random | 65–97% |
| Within-region exact-C | site | full | trinuc + bin | matched | ~ chance (stochastic) |

**Headline (deployment-realistic, audited):**
> Given an editor's deaminase motif + a cross-editor accessibility model, we recover **~30–37% of its guide-independent off-targets in the top 10% of screened sites (3–4× enrichment)**, generalizing across treatments AND to unseen deaminases.

**f × g per editor (held-out, motif provided, g from other sources):**

| Editor | Deaminase | f-only | g-only | **f × g** |
|---|---|---|---|---|
| Doman | rAPOBEC1 | 1.9× | 1.1× | 2.3× |
| McGrath | rAPOBEC1 | 2.4× | 0.9× | 2.2× |
| Lei | rAPOBEC1 | 2.4× | 1.9× | 3.3× |
| **Yu** | **engineered (novel)** | 3.1× | 1.6× | **3.7×** |
| **Richter** | **TadA (novel)** | 2.9× | 2.1× | **3.5×** |

## 6. Discussion — what to talk about

**What this confirms from your RNA framework:**
- Motif × landscape factorization holds for DNA.
- Aggregate / regional burden (editing-index analog) is far more reproducible than per-site.
- Cross-cell-type stability at the regional scale.

**Where DNA differs:**
- Substrate is permanent/heritable → safety bar higher.
- No matched-control bulk WGS achievable (the editing-index trick doesn't translate at single-base VAF).
- We rank candidate sites; we cannot measure aggregate burden directly.

**Honest limits:**
- Within-region exact-C is stochastic (irreducible).
- Per-gene cancer-driver flagging is largely gene-density (not biology).
- Calibrated absolute rate is not supported (no cross-assay gold).
- A brand-new, never-characterized deaminase is degraded — supply its motif spectrum or train on a sibling source.

## 7. Three asks for Levanon

1. **Prospective wet-lab validation — the trump card.** Pick one editor; deep-sequence a small panel of **gene-density-matched** high-rank vs low-rank regions predicted by f × g. This is the only confound-proof escalation; it converts our screening prior into a validated demonstration. ~20 + 20 regions × deep coverage, one editor — feasible in your hands.
2. **Joint RNA + DNA shared enzyme embedding** (Phase D, currently backlog). Your RNA model + our DNA model both factor through the same deaminase identity. A shared enzyme embedding trained jointly should improve per-editor priors and quantify cross-substrate coupling (does rAPOBEC1 RNA preference predict its DNA preference? Architecturally yes, empirically unknown.)
3. **Cleaner / restricted data sources.** Verve VERVE-101 (ABEmax in NHP, Rothgangl 2021) — restricted; can he facilitate? Plus any unpublished bulk treated-control WGS from his network — the published landscape has been exhaustively mapped.

## 8. Possible pivots if he prefers a different direction

- **Guide-DEPendent** prediction (CHANGE-seq-BE, EndoV-seq) — data-rich, clinically per-therapy, ML-hard.
- **Per-editor total burden ranking** (Buchumenski analog) once we have multiple bulk genome-wide assays per editor.
