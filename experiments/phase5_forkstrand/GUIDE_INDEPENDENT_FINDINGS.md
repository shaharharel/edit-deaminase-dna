# Guide-Independent Base-Editor DNA Off-Targets — Findings & Model-Card

**Status:** investigation COMPLETE — negative, airtight (multiple independent QA passes + a region-matched closure control).
**Date:** 2026-07-19 (24h autonomous run). **Scope:** guide-*independent* (deaminase-intrinsic) DNA off-targets only.
**Datasets:** GOTI (Zuo 2019, mouse 2-cell embryo WGS, ~30×) · Doman (2020, human 24-clone WGS, ~30×) · Lei (2021, human 293T Detect-seq).

---

## Bottom line (one paragraph)

Guide-independent base-editor (rAPOBEC1 CBE) DNA off-targets **cannot be isolated from the cell's endogenous
APOBEC3 background, nor predicted at the single-site level, in 30× WGS.** The one real, reproducible biological
signal — editing prefers the replication-fork-exposed (lagging-strand-template) ssDNA at TpC motifs — is the
**known APOBEC lagging-strand rule and is SHARED**: the cell's own APOBEC3 obeys it identically, so it does not
separate editor sites from background and is not editor-specific. As a per-site *predictor* it is ~chance
(AUC ≈ 0.50–0.52, both species). It survives only as a mechanism-grounded **risk covariate** (one input to a
genome-wide risk map), **not** a standalone safety gate. The root cause of the negative is a hard sensitivity floor:
true diffuse low-VAF guide-independent editing sits **below the detection floor at 30×**, so the called-variant set
is sequencing artifacts + shared background, and the editor's contribution is a small dispersed increment no
statistic can resolve. Revisiting requires **100×+ error-corrected (duplex/NanoSeq) sequencing** — new data, not
new analysis.

---

## Two questions asked, both answered NEGATIVE

### Q1 — Editor-vs-endogenous separability (can we attribute a site to the editor?)  → NO
Tested on five orthogonal axes; all negative or uninterpretable. The decisive control holds Cas9 constant to isolate
the deaminase (Doman BE4 = deaminase+nCas9 vs nCas9-only = Cas9, no deaminase), at clone level.

| Axis | Result | Verdict |
|---|---|---|
| C→T burden / editing index | editor only **+4%** over nCas9; motif-degenerate | not separable |
| Replication-strand coupling (clone-level BE4 vs nCas9) | Mann-Whitney **p=0.29** (p=0.54 at n=16 adding YE1-BE4); the *weaker* YE1-BE4 editor showed **less** coupling (0.85) than deaminase-free nCas9 (1.39) — wrong direction | shared, not editor-specific |
| C→G / UGI signature | called C→G (3425) **exceeds** called C→T (2891) at Parent-clean TpC — physically impossible for a deaminase → the ≥3-read call set is **sequencing artifact** (C→G≈C→A≈C→T floor) | dead |
| Pentamer YTCA/RTCA | BE4 1.060 vs nCas9 1.030 (~3%), confounded by differential APOBEC3 load | inconclusive/confounded |
| Clustering / kataegis | ~360 called C→T/clone → expected same-strand pairs within 1 kb **<<1/clone** | power-dead (NO-GO) |

**Why they all fail (one root cause):** the true editor off-target signal is diffuse and low-VAF, **below the ≥3-read
detection floor at 30×**. What is measurable is artifacts + the shared endogenous-APOBEC3 background; the editor's
increment is a small, dispersed addition unresolvable by any summary statistic on ~360 events/clone.

### Q2 — Per-site predictability (can we build a location risk-map?)  → ~CHANCE
Feasibility = predict edited vs *detectability-matched* unedited sites, leave-one-chromosome-out, vs a callability
baseline, permutation-null.

| Dataset | Matching | Baseline | +strand | Verdict |
|---|---|---|---|---|
| GOTI (mouse) | NN on trinuc/GC/TSS/rep (region-aware) | 0.43 (≈chance after matching) | 0.524, macro-Δ **+0.093**, perm p=0.005 | real strand lift, but absolute **~0.52** = weak |
| Doman (human) | coverage-matched | 0.497 | magnitude Δ≈0 (p=0.52); +on_exposed 0.520 (+0.023) | ~chance; +0.023 was region confound → |
| Doman (human) | **region-matched** (coverage × RFD-sign × \|RFD\|-tertile) | 0.510 | +on_exposed 0.504 (**no lift over baseline**) | **+0.023 collapses → airtight chance** |

The GOTI strand lift is real and detectability-independent, but absolute AUC tops out at ~0.52 — a diffuse-process
ceiling, not a usable site predictor. The one human positive number (+0.023) was replication-domain confound and
vanished under region-matched negatives.

---

## What IS real (and worth keeping)

The **replication-fork lagging-strand rule** for guide-independent editing, measured genome-wide via OK-seq RFD:
- callability-immune C>T-vs-G>A × sign(RFD) split; composition-corrected (ascertainment-matched null).
- **GOTI mouse ~2.47** · **Doman human 1.41 pooled / ~2.1 at high |RFD|** (genic universe composition-flat) ·
  **Lei human 1.32 pooled / 1.71 at high |RFD|** (whole-genome universe, composition-corrected).
- Dose-responds with fork polarization (|RFD| tertiles). This is the `f(sequence:TpC) × g(structure:fork-exposure)`
  DNA analog of the RNA model's motif×structure — it works mechanistically but is **shared and weak**.

**Use:** a mechanism-grounded, transferable **risk covariate** (shared across CBE constructs and endogenous APOBEC3 —
a generalization advantage for a safety feature), fed into a larger DeaminaFormer-DNA model. **NOT** an
editor-specific signal, **NOT** a per-site predictor, **NOT** a standalone gate.

---

## Confounds policed (and how)
- **Callability/detectability:** feasibility used detectability-matched negatives (not trinuc-only) + a callability
  baseline the strand model had to beat; the human positive collapsed under stricter region-matching.
- **Composition:** ascertainment-matched (not genome-wide) reference-composition null — negligible for genic Doman
  (~1.02 flat), real for whole-genome Lei (corrected out). Genome-wide null over-corrected (a caught error).
- **Clonal-APOBEC3 load:** the C>T-vs-G>A×RFD split is symmetric to editing load; clone-level tests avoid pseudoreplication.
- **Between-culture (Lei NT vs Empty):** replicated (rep2) but retired as editor-specific by the Cas9-constant
  clone-level control (NT−Empty confounds deaminase with Cas9/vector).
- **Artifact floor:** the C→G channel positively *demonstrated* the ≥3-read call set is artifact-dominated.

## Independent replication
The sister RNA project independently concluded human guide-independent WGS burden is a valid **index** but **not
per-site concentrable** (motif-degenerate, strand-flat, APOBEC3-co-reproducible) — reached via a *different* assay and
*different* failure modes, so this is genuine convergence, not groupthink.

---

## Forward paths (all require a directive or new data)
1. **100×+ duplex / NanoSeq error-corrected sequencing** — the only way to push the detection floor below single-read
   errors and recover the faint diffuse editor edits. New data; would re-open Q1.
2. **Pivot to guide-DEPENDENT off-targets** — Cas9-binding-driven; prior work found real, concentrable seed/PAM signal
   that beats the mismatch-count baseline. The genuinely predictable DNA direction, with data in hand.
3. **Integrate the fork-strand covariate** into DeaminaFormer-DNA as one interpretable replication-topology feature
   (cheap; does not by itself justify Evo+chromatin).

## Reproducibility
Scripts: `experiments/phase5_forkstrand/` (laptop `/tmp/poc_dna/`, box `/data/scripts/`). Key: `clone_level_be4_vs_ncas9.py`
(decisive editor-specificity control), `goti_feasibility_hardened.py` + `doman_feasibility_regionmatched.py`
(feasibility + closure), `lei_composition_null.py` / `lei_covered_null.sh` (composition nulls), `cg_channel_test.py`
(artifact-floor). Full audit trail: `iterlog_local.md` (mirrored `gs://ai-temp/apobec-genome-cache/poc24_results/`).
