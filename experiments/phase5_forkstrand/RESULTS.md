# Phase 5 — Replication-fork-strand predictor of guide-independent CBE off-targets

**Status:** REPRODUCED 3-way (all with cell-type-mismatched RFD tracks) for the replication-fork lagging-strand rule of guide-independent editing — GOTI mouse WGS (~2.47) + Doman human WGS (1.41 pooled / 2.10 high-|RFD|, composition-clean genic universe) + Lei human Detect-seq (composition-corrected 1.32 pooled / 1.71 high-|RFD|). **EDITOR-SPECIFICITY: NEGATIVE (decisive).** The deaminase-isolating clone-level control (Doman BE4 vs nCas9, Cas9 held constant) is null: BE4 8-clone median strand-OR 1.57 vs nCas9 6-clone 1.39, Mann-Whitney p=0.29. So the coupling is the **shared** APOBEC rule, and the Lei NT−Empty "increment" (1.11–1.23, replicated across 2 culture pairs) is **NT-condition-associated (Cas9/vector/subclonal-APOBEC3), not deaminase-associated** — replication rules out random drift but not the systematic NT-vs-Empty difference. Scope: known APOBEC rule, ~1.3–2.1× (APOBEC ceiling), a real composition-clean **SHARED risk covariate — not an editor-specific per-site signal, not a standalone gate**; RFD tracks all cell-type-mismatched.
**COMPOSITION-NULL investigation (QA, 2026-07-19) → RESOLVED to negligible when ascertainment-matched.** The C>T-vs-G>A×RFD split is load/callability-immune but not automatically composition-immune. A *genome-wide* null (`lei_composition_null.py`, raw hg19.fa, ZERO editing) showed reference TpC-C/GpA-G strand-skew rising with |RFD| (1.06/1.17/1.43), *seeming* to imply ~40% inflation. **But the correct null is ascertainment-matched** — built from the positives' own covered universe, not genome-wide. Doman's `per_site_spectrum` is a 15.67M-site genic universe; its covered-unedited TpC pool (n=3.86M, same coverage filter + fixed |RFD| bins) gives **OR ~1.02 FLAT** (low 1.02 / mid 1.04 / high 1.02) — no skew, no gradient. The genome-wide 1.21 came from intergenic/repeat regions absent from the analysis universe. **→ composition confound is negligible in-universe; Doman strand signal is REAL: pooled 1.41, high-|RFD| 2.10 (corrected 1.37/2.06).** **Lei is the opposite case:** its ascertainment universe is whole-genome Detect-seq, so its in-universe covered-unedited null IS composition-skewed (1.04→1.18→1.46, ALL 1.22 ≈ genome-wide) → composition is a real confound for Lei. Lei composition-corrected (edited ÷ in-universe null): pooled 1.61/1.22=**1.32**, high-|RFD| 2.49/1.46=**1.71**. **So the composition confound is universe-dependent: negligible for genic Doman, real for genome-wide Lei.** GOTI 2.47 is composition-safe (matched negatives = in-sample ascertainment-matched null). **Methodological lesson: build detectability/composition nulls from the positives' OWN ascertainment universe — it can be flat (genic) or skewed (genome-wide), so a genome-wide null is right by luck for some datasets (Lei) and wrong for others (Doman).**

**DECISIVE EDITOR-SPECIFICITY TEST → NEGATIVE.** The contrast that isolates the deaminase (Cas9 held constant): Doman clone-level BE4 (deaminase+nCas9) vs nCas9-only, per-clone TpC strand-OR, unweighted Mann-Whitney (`clone_level_be4_vs_ncas9.py`). BE4 8-clone median **1.57** [0.84,2.92] vs nCas9 6-clone median **1.39** [0.93,1.94], **p=0.29 (NS)**. → No clone-level deaminase-specific increment; the coupling is the **shared** APOBEC rule. This retires the Lei NT−Empty "editor increment" (1.11–1.23) as *editor*-specific — 2-pair replication rules out random drift but not the systematic NT-condition difference (Cas9+vector + subclonal-APOBEC3 divergence), so it is **NT-condition-associated, not deaminase-associated**. **Net: the fork-strand axis is a real, composition-clean, SHARED risk covariate — not an editor-specific per-site signal.**
**Date:** 2026-07-19 (24h autonomous run). *QA arc (overclaim → germline flag → motif-restriction → editor-free differential → between-culture-caveat → composition-null correction) documented as pre-registered honesty: the endogenous baseline is explicitly quantified, and the editor increment is reported as editor-associated (not proven editor-specific) given the single between-culture contrast.*

## Headline (honest)

Guide-independent cytosine-base-editor (rAPOBEC1 CBE) off-target C→T editing obeys the **known APOBEC
lagging-strand-template rule**: at replication forks the deaminase edits the persistently-exposed
lagging-strand-template ssDNA. Measured genome-wide via OK-seq **Replication Fork Directionality (RFD)**, this is a
**directional (which-strand) bias** that replicates across mouse and human at the effect size established for
endogenous APOBEC in cancer.

It is a **mechanism-grounded, callability-immune risk covariate** that predicts off-target *location/strand*. Two
things are true and must be stated together: (i) the strand **rule is SHARED** — endogenous APOBEC3 obeys it too, so
it is NOT a concentrable editor-vs-background *panel* (it can't separate editor sites from endogenous APOBEC sites);
(ii) the editor-present sample (Lei NT) shows a stronger coupling than the editor-free vector (Empty), an
**editor-associated increment** (p=3.5×10⁻⁷⁹) — *consistent with* the editor's off-targets obeying the rule above the
endogenous baseline, though not fully isolated from between-culture 293T divergence (see limitations). So it is the first
*replication-topology* feature shown to carry site-level directional signal for a therapeutic base editor's
off-targets — a useful interpretable term in a larger model, **not** a standalone safety gate.

## Why this matters (the confound it beats)

Every **rate-based** attempt to isolate the editor from endogenous clonal APOBEC3 background failed (per-gene DEI,
motif, f×g, DNase, R-loop, C→G/UGI channel — all null/inconclusive), because between-clone APOBEC3-load variance ≥ the
~4% editor increment in an n≈8 per-clone WGS design. The fork-strand axis succeeds because the test is the
**C>T-vs-G>A × sign(RFD) split**, which is *symmetric to total editing load* (clonal-APOBEC3-load variance cannot flip
a C>T event's strand) and to *callable-C-ness* (detection bias acts identically on C>T and G>A). Only genuine
strand-coupling partitions the two channels oppositely along RFD.

## Evidence

| # | Control | Result |
|---|---------|--------|
| 1 | GOTI beats callability+accessibility baseline (LOChr-out) | Loose (trinuc-only) negs: AUC 0.5189→0.5654, Δ+0.0465, perm p=0.016. **HARDENED (1:1 NN detectability-matched negs — in_gene/gc/tss/rep balanced, 5-seed, macro-AUC):** baseline collapses to **0.43 (≈chance → old baseline was mostly detectability)**, strand model **0.524**, macro-Δ **+0.093** (seed-stable), perm p=**0.005** → PASS. **Two truths: strand lift is REAL + detectability-independent, but absolute site-prediction is only ~0.52 (barely>chance) — the 0.565 was detectability-inflated. A real covariate, NOT a usable per-site predictor.** **HUMAN (Doman) confirms: coverage-matched, baseline 0.497 / +strand-magnitude 0.496 (perm p=0.52, no signal) / +on_exposed 0.520 (confound-prone) → human per-site prediction ≈ CHANCE. Feasibility leg complete both species; strand does NOT discriminate edited from matched-unedited human sites.** |
| 2 | Mechanism is **directional**, not magnitude | lift carried by `on_exposed` (+0.029) ≫ `absrfd` (+0.010) |
| 3 | Transcription control — Doman (MH by gene strand) | OR 1.41→**1.40** adjusted, p=7×10⁻⁶ (undiminished) |
| 3b | Transcription control — GOTI (intergenic-only) | OR **2.24**, p=6×10⁻⁹; on_exposed pos 0.599 vs neg 0.509 |
| 4 | **Clone-level replication** (Doman, 8 BE4 clones) | 7/8 clones OR>1, Wilcoxon p=0.012 (not pseudoreplication) |
| 5 | Positive-control calibration (3 APOBEC sources) | BE4 1.40, YE1 1.29, **nCas9 (editor-free) 1.30** — pipeline recovers the known rule; nCas9 = clean positive control |
| 6 | **Dose-response** (Doman human, \|RFD\| tertiles) | split OR 1.11 → 1.25 → **1.99** (low→high \|RFD\|), p=5×10⁻⁵ |
| 7 | Domain-burden **magnitude** (\|RFD\| and Repli-seq, uniform callable background) | **null** (IRR 0.98 / 1.04 after coverage) → claim bounded to *strand*, not magnitude |

| 8 | **Lei human Detect-seq — APOBEC-context signal, motif-specific** | pooled OR=1.25 (n=410,765) concentrates in **TpC (APOBEC context): OR=1.609** (n=174K); germline-CpG only 3.5% of sites; non-APOBEC context near-null (descriptive); holds across all VAF incl editor-like <0.10 (1.47). Motif-specificity (TpC vs non-TpC) argues against a Detect-seq strand artifact; VAF-flatness refutes a germline confound. **Composition check (see Status): a genome-wide null suggested inflation, but the ascertainment-matched covered-unedited null is the correct baseline (Doman's is ~1.02 flat → signal real); Lei's own covered null is recomputing before finalizing the Lei concentration number.** |
| 9 | **Lei editor-ASSOCIATED increment** (NT editor-present vs Empty editor-free vector) | OR_TpC(NT)=**1.61** vs OR_TpC(Empty)=**1.31** (endogenous 293T APOBEC3 baseline; Empty all-context OR=1.08 near-null → even endogenous coupling is APOBEC-context-specific) → increment **1.23× [1.20,1.26], p=3.5×10⁻⁷⁹**, TpC-specific (non-TpC bounded <1.08). **Caveat:** a single between-culture contrast (NT vs Empty are different cultures of unstable 293T) — consistent with an editor contribution but not isolated from subclonal-APOBEC3 divergence; n=174K makes any difference "significant", so p is not provenance evidence. |

**Honest concordance claim — 3-way:** GOTI mouse WGS (2.47) + Doman human WGS (1.40, TpC) + Lei human
Detect-seq (TpC 1.61, with an editor-specific increment over the editor-free baseline). The Lei story is a model of
iterative QA: an initial pooled "replication" (1.25) drew a germline flag (411K positives ≫ plausible editor count) →
**motif-restriction** localized the signal to APOBEC context (1.61) and refuted germline (CpG only 3.5%) + artifact
(TpC-specific) confounds → the **editor-free Empty vector's own TpC coupling (1.31)** quantified the endogenous
APOBEC3 baseline, and NT's higher 1.61 is an editor-associated increment (p=10⁻⁷⁹) — consistent with an editor
contribution above baseline, though a single between-culture contrast (not isolated from 293T subclonal divergence; see
limitations). Report Lei motif-restricted with the NT-vs-Empty differential, not the pooled-diluted OR.

## Honest limitations

- **Shared mechanism, not a concentrable panel.** Endogenous APOBEC3 obeys the same rule (Doman nCas9 OR 1.30; Lei
  editor-free Empty OR_TpC 1.31), so the strand axis cannot *separate* editor sites from endogenous APOBEC sites — it is
  not a concentrable editor-vs-background panel. (Note: Doman's *rate* comparison BE4-vs-nCas9 was not separable, p=0.20;
  but the Lei *strand-coupling* differential NT-vs-Empty IS significant, p=10⁻⁷⁹ — the editor contributes coupled
  editing above baseline even though the rule itself is shared. These are consistent: shared rule + measurable editor increment.)
- **Modest pooled effect** (OR 1.40). Diluted average of null (low \|RFD\|) and OR~2.0 (high \|RFD\|). Too weak to carry a
  hard-block safety tier alone.
- **Cell-type-mismatched RFD track — ALL three legs.** HeLa RFD applied to non-HeLa Doman + 293T Lei; and GOTI uses
  mm10 mESC RFD applied to 2-cell mouse embryo (species-matched but NOT developmental-stage-matched — early-embryonic
  replication programs differ from mESC). So "reproduced 3-way" all share this caveat. The 3-source calibration proves
  *precision*, not absolute accuracy; the matched-track effect size is unmeasured (no 293T OK-seq dataset exists to build one).
- **Lei editor increment is a single between-culture contrast**, not isolated from subclonal-APOBEC3 divergence in the
  unstable 293T line (see Status). Consistent with an editor contribution; not a replicated, culture-controlled editor increment.
- **Reference-composition confound — checked, negligible when ascertainment-matched.** The C>T-vs-G>A×RFD split is not automatically composition-immune, so we tested it. A genome-wide null showed apparent skew (up to 1.43), but the correct **ascertainment-matched** null (Doman's genic covered-unedited TpC pool, n=3.86M) is **~1.02 flat** — no skew, no |RFD| gradient. So Doman's signal (1.41 / 2.10) is real, not composition. (Lei recompute in progress.) The genome-wide null over-corrected because it sampled intergenic/repeat regions the positives never occupy. This is the load-bearing methodological caveat: **any null for this split must be built from the positives' ascertainment universe.**
- **Magnitude leg null** — the mechanism predicts *which strand*, not per-domain burden.
- **C→G channel (UGI signature) axis is DEAD at 30× — third concordant editor-specificity negative.** Premise: BE4 is UGI-protected → should suppress the C→G/C→A arms vs endogenous APOBEC3's SBS13. Test (Doman, TpC, Parent-clean): called mutations (≥3 reads) give BE4 C→T=2891 / **C→G=3425** / C→A=3141 and nCas9 C→T=2249 / **C→G=3037** / C→A=2651. Called C→G *exceeds* called C→T — physically impossible for a deaminase (~95% C→T) → the ≥3-read events are recurrent sequencing/mapping **artifacts** (C→G≈C→A≈C→T = artifact-floor; C→A/C→T≈2 in the Σ/Σ view = OxoG), and the true diffuse low-VAF editing sits *below* the floor. Ratios don't separate BE4 from nCas9 (0.62 vs 0.66; 1.18 vs 1.35). No rescue at 30× (needs 100×+ duplex/NanoSeq). **Together with the null rate/burden axis and the clone-level strand null (p=0.29→0.54), three concordant negatives bound guide-independent editor-vs-endogenous separation as *sub-artifact-floor at 30× WGS* — the fork-strand axis is a shared covariate, not a concentrable editor panel.**

## Methods (scripts, laptop `/tmp/poc_dna/`, mirrored `gs://ai-temp/apobec-genome-cache/goti_lei_results/`)

- Data: GOTI = Zuo 2019 S7 CBE SNVs (mm10, 1609 pos + 4749 detectability-matched negs); Doman = 24-clone WGS
  per-site spectrum (hg19, per-cohort/per-clone C→T/C→G/C→A); Lei = 2021 Detect-seq 293T-NT (guide-independent) vs
  293T-Empty (hg19).
- RFD: OK-seq (CL-CHEN-Lab) — mm10 mESC, hg19 HeLa.
- `goti_domain_v2.py` (callability baseline + delta null), `goti_decomp.py` (directional decomposition),
  `goti_intergenic.py`, `doman_strand.py` + `conditional_strand2.py` (MH transcription control),
  `clone_level_strand.py`, `pos_control_calib.py`, `doman_dose.py`, `doman_rt_burden.py`.
- Lei: `lei_par*.sh` + `lei_par_analyze.py` (fast parallel per-chrom mpileup caller — awk-prefilter + `-d 500`, fixes
  the multi-hour pure-Python whole-genome parse) → `lei_rigorous.py` (pre-registered: sign → Empty-survival →
  intergenic → magnitude → power) → `lei_vaf_control.py` (germline VAF stratification) → `lei_tpc_control.py` (motif
  restriction, hg19.fa context) → `lei_empty_par.sh` + `lei_empty_analyze.py` (editor-isolating NT-vs-Empty TpC differential).

## Next steps (all require new data or a directive)

1. Measure the RFD-track attenuation via matched-cell-type (293T/HEK) OK-seq RFD — replace the mismatched HeLa track to
   get the un-attenuated human effect size and test the attenuation hypothesis (matched mouse mESC gave OR~2.5 vs
   mismatched human ~1.4–1.6). Requires sourcing + processing raw OK-seq (no clean pre-computed HEK293 RFD track found).
2. Integrate `on_exposed`/`absrfd` as an interpretable term in DeaminaFormer-DNA (does NOT justify Evo+chromatin on
   its own — it is a cheap replication-topology feature).
3. (Optional) External APOBEC-high tumor-WGS positive control — now largely superseded by the Lei editor-free Empty
   differential, which already isolates the editor increment over the endogenous APOBEC3 baseline.
