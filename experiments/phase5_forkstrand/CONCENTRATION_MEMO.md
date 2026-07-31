# What f(sequence) × g(structure) can and can't do for guide-independent CBE liability

**Question (from PI):** can a per-site risk model — `f(sequence) × g(structure)` — rank the diffuse guide-independent
off-target editing into a shorter set of "highly-probable" sites, with 2–4× enrichment of editing index in the top-N?

**Short answer:** the model *shape* is right and every component is *real and editor-detectable*, but each is *weak*,
and no combination reaches a 2–4× enriched panel — because the guide-independent editor signal is a small excess over a
large, same-motif background (endogenous APOBEC3 + sequencing noise). It is a **burden index + directional risk
covariate**, not a concentrable panel.

## The components we tested (all clone-level / LOChr-out / matched-control where applicable)
| Component | Type | Real & editor-specific? | Effect size |
|---|---|---|---|
| Motif TpC | f(seq) — substrate | yes, shared with APOBEC3 | ~1.16× (weak; "motif-degenerate") |
| Pentamer YTCA/RTCA | f(seq) — deaminase fingerprint | **yes, editor-distinct** (BE4 1.060 vs nCas9 1.030, clone-level p<0.001) — but not fully isolated from differential A3A load | ~3% (tiny) |
| Fork-strand (RFD on_exposed) | g(structure) — ssDNA exposure | yes (direction), shared | ~1.4× on *strand*; ~1.0× on *site-rank* (human) |

## What the metrics showed
- **Site-classification** (vs detectability-matched negatives, GOTI): fork-strand delta +0.047 AUC (beats callability
  baseline, bootstrap p<1e-4 + perm p=0.016) → real but modest; top-10% captures 14% of editor sites (lift 1.4×).
- **Editing-index enrichment** (top-N vs genome, Doman): **~1× (no enrichment)** — the C>T *rate* is dominated by
  background (germline C>T-at-CpG + sequencing error); the editor's ~1.16× TpC excess is buried, so ranking can't
  concentrate it.
- **Fork-strand isolated from coverage** (Doman): ~1.0× site-concentration — the apparent ranking was mostly
  detectability (editor sites need coverage to be *called* → higher-coverage sites rank higher).

## Why this is structural, not a modeling gap
The strand *ratio* test works (OR 1.40) because it is immune to the background *rate* — it asks "of the C>T that
happened, which strand?" But ranking/enrichment asks about the *amount*, which is background-dominated. Direction
survives the noise; magnitude does not. Stacking three weak signals (1.16× × 1.03 × 1.4×-on-strand) does not amplify a
signal that is <20% of what is measured.

## What the model DOES deliver (honest, useful)
1. A per-sample **guide-independent burden index** (total deaminase-intrinsic liability).
2. A mechanism-grounded **directional risk covariate** — motif × pentamer-fingerprint × fork-strand — that reliably
   predicts *which strand / which context* is at higher risk. An interpretable term for a larger model (the
   DeaminaFormer f(seq)×g(structure) shape), with modest effect size = the true biology.

## What it does NOT deliver
A 2–4× enriched, short "highly-probable mutations" panel from guide-independent editing. Concentration ceiling ≈ 1.4×
(site-density) / ~1× (editing-index). The process is diffuse and motif-degenerate with the cell's own APOBEC3.

## Strategic implication
- **Guide-independent (deaminase-intrinsic):** report as an **index + directional covariate**. The pentamer fingerprint
  confirms rAPOBEC1 is *distinct* from endogenous APOBEC3 — scientifically interesting, but too subtle to concentrate.
- **A short high-yield panel** is achievable only for **guide-DEPENDENT** off-targets (gRNA/protospacer-homology driven,
  ~4–5× validation-panel reduction shown earlier) — but that is a property of the gRNA targeting, NOT the deaminase.
- To push guide-independent toward a panel would require a fundamentally cleaner editor-specific readout (deep
  enrichment + site-by-site matched no-editor subtraction), not more features on the current data.

---
## ADDENDUM (assay matters): Detect-seq concentrates the *directional* signal better than WGS
The "~1× enrichment / no panel" conclusion above was measured on **WGS (Doman)**, where editor sites are buried in
background. On **Lei Detect-seq** — the enrichment assay that pulls down edited (deoxyU) sites — the directional
strand-coupling concentrates with `TpC (motif) × high-|RFD| (fork polarization)`:

- Raw fork-strand OR by |RFD| tertile *within TpC* (411K called sites): low 1.19 → mid 1.45 → high 2.49.
- **Caveat:** this concentrates the *directional/strand-coupled* editing (which exposed-strand APOBEC-context site),
  NOT the raw editing *amount* (VAF flat ~0.27 across strata); and editor-vs-endogenous is not separated by it.

## ADDENDUM-2 (QA composition-null correction, 2026-07-19): the raw 2.46× was ~40% reference composition — "2–4×" REFUTED
A scientific-analyst QA flagged that the C>T-vs-G>A×sign(RFD) split, while load- and callability-immune, is **NOT
reference-composition-immune**: replication-initiation zones (high |RFD|) are strand-composition-skewed, so the TpC-C vs
GpA-G *reference denominators* are asymmetric w.r.t. RFD sign. A pure-genome null (`lei_composition_null.py`, hg19.fa,
**zero editing**, 4.85M sites) reproduces the gradient with no editing at all: composition-OR low 1.06 / mid 1.17 /
high 1.43 / ALL 1.21.

Composition-corrected (edited ÷ composition, fixed |RFD| cuts):

| |RFD| | raw edited | composition (no editing) | **corrected** |
|---|---|---|---|
| low | 1.19 | 1.06 | **1.12** |
| mid | 1.45 | 1.17 | **1.25** |
| high | 2.49 | 1.43 | **1.75** |
| ALL TpC | 1.61 | 1.21 | **1.34** |

> ⚠️ **SUPERSEDED — the genome-wide null above was the WRONG baseline (see ADDENDUM-3).** The `1.34/1.75`
> "corrected" numbers and the "2–4× withdrawn" conclusion in this block were computed by dividing by a *genome-wide*
> composition null, which sampled intergenic/repeat regions the edited positives never occupy. The correct null is
> **ascertainment-matched** (built from the positives' own covered universe), and it changes the answer. Read
> ADDENDUM-3.

## ADDENDUM-3 (correction to ADDENDUM-2, 2026-07-19): the composition confound is NEGLIGIBLE when the null is ascertainment-matched
The composition null must be built from the population the positives are **drawn from**, not genome-wide. Test case
(Doman, where the full covered universe is available locally): Doman's `per_site_spectrum` is a 15.67M-site **genic**
universe; its **covered-but-unedited** TpC-C vs GpA-G pool (n=3.86M, same coverage filter, same fixed |RFD| bins as the
edited analysis) gives strand **OR ≈ 1.02, FLAT across all |RFD| tertiles** (low 1.02 / mid 1.04 / high 1.02) — no skew,
no gradient. The genome-wide 1.06→1.43 gradient came from regions **absent** from the analysis universe. So for Doman
the composition confound is negligible and the **strand signal is real** (pooled 1.41, high-|RFD| 2.10; corrected
1.37/2.06). **Lei is being recomputed** with its own Detect-seq covered-unedited null (`lei_covered_null.sh`); the
ADDENDUM-2 "corrected 1.34/1.75" is withdrawn pending that in-universe null.

**The NT−Empty "increment" is NOT editor-specific (retired by the decisive clone-level control).** The Lei increment
(rep1 1.23; rep2 1.11 [1.08,1.15]) replicates across two culture pairs — but replication only rules out *random*
drift, not the *systematic* NT-vs-Empty difference (NT carries a Cas9 + sgRNA vector Empty lacks, plus possible
subclonal-APOBEC3 divergence — both replicate across pairs). The contrast that actually isolates the deaminase (Cas9
held constant) is Doman clone-level **BE4 vs nCas9-only**: BE4 8-clone median strand-OR **1.57** vs nCas9 6-clone
**1.39**, Mann-Whitney **p=0.29 (NS)**. → **No clone-level deaminase-specific increment.** The coupling is the shared
APOBEC rule; the Lei increment is NT-condition-associated, not deaminase-associated.

**Net on the panel question (final):** the directional strand-coupling is real and composition-clean, but (i) modest
(composition-corrected ~1.3–2.1× at high |RFD|) and (ii) **SHARED, not editor-specific** — so it cannot separate
editor sites from endogenous-APOBEC3 sites and does not yield a per-site editor panel. Lei composition-corrected
(in-universe null 1.22/1.46): pooled **1.32**, high-|RFD| **1.71** — below the "2–4×" target, and not editor-specific
anyway. Magnitude enrichment (VAF) and editor-vs-endogenous *rate* separation remain out of reach. GOTI's 2.47 is
composition-safe (matched negatives). **Conclusion: fork-strand is a real shared risk covariate, not a concentrable editor panel.**
