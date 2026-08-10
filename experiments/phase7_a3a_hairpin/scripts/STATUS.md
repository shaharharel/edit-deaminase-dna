# A3A / PCAWG run — STATUS
_updated 2026-08-09 23:20 by autonomous monitor_

## Hypothesis under test
Train a sequence/structure classifier on ENDOGENOUS APOBEC3A mutations from cancer
genomes (PCAWG), then ask whether it predicts (a) endogenous A3A activity in
base-edited clones and (b) A3A-class base-editor off-target burden.
Governing idea: hairpins make ssDNA exposure SEQUENCE-DETERMINED, which is the one
mechanism that can transfer; cell-state ssDNA proxies have all been null.

## DONE
- S1.1 PCAWG open MAF (925MB, no DACO) downloaded + gzip-verified
- S1.1 parsed 21,606,842 usable SNVs / 1,830 donors, pyrimidine-oriented, hg19
- S1.1 per-donor APOBEC enrichment (1,776 donors) + -2 base for YTCA/RTCA
- S1.1 donor ranking: 1,637 donors >=50 TCW, 3,218,946 TCW mutations total
      A3A-like (tcw_frac>=0.20 & ytca_frac>=0.55): 94 donors, 348,216 TCW muts

## RUNNING
- a3a-s1-stage2 : 100k positives + 10x trinuc-matched/donor-matched negatives,
                  hairpin scan (stem 3-10, loop 3-9), local GC -> feat/a3a_trainset.npz
- a3a-s2-w1/w2/w3 : PRJNA1042830 clonal WGS download+align, 3 workers x 9 threads
                  w1 Parent, Y130G-clone1, VA-clone1
                  w2 nCas9-clone1, Y130G-clone2, VA-clone2
                  w3 nCas9-clone2, YE1-clone1, YE1-clone2

## KNOWN CAVEAT (important)
PRJNA1042830 = Nat Chem Biol 2024 PMID 38553609, haA3A-CBE variants ENGINEERED
for near-background off-target activity. So the Y130G/VA editor arms may carry
little editor-specific burden -> low power for the editor-burden test.
PRE-REGISTERED: we expect NO enrichment on haA3A editor-specific calls. This is
recorded BEFORE looking, so neither outcome can be retrofitted as success.
The endogenous-A3A test (Parent + nCas9 clones) is unaffected and remains primary.
Wild-type A3A-BE3/A3A-BE4max clonal WGS would be the high-burden arm — still being
hunted by the a3a-editor-survey agent.

## NEXT
1. Train CatBoost/GB on hairpin features, held-out CHROMOSOME. Baseline to beat.
2. DNA LM (HyenaDNA / Caduceus / NT-500M+LoRA) on V100 ai-gpu. NOT Evo-7B.
3. Score Doman clones + PRJNA1042830 clones; nCas9 calibrator, 3x bar.

## DISCIPLINE (non-negotiable)
- Negatives trinuc-matched (unmatched -> meaningless AUROC ~0.68 = motif-learning)
- RANDOM baseline printed beside every enrichment number
- Held-out by CHROMOSOME, never random split
- nCas9/Parent calibrate every editor claim; editor increment must be >=3x

## RESULT — S3 classifier v2 (2026-08-09 23:45), CORRECTED build
BUG FOUND+FIXED: stage2 sampled negatives plus-strand-only while positives are
pyrimidine-oriented (50,069/100,000 minus-strand). revcomp(TCT)=AGA has flanks
A..A vs TCT's T..T => sequence model separated minus-strand TCT positives
perfectly. Symptom: top-1% enrichment exactly 11.000x = the ceiling, top slice
100% TCT. v1 numbers (AUROC 0.735, 11x) are DISCARDED.
stage2b fix: negatives matched on (trinuc, STRAND), all features computed in
strand-oriented windows. Leak checks now pass (focal base C 1.0000/1.0000;
offset -1 and +1 ACGT distributions identical pos vs neg).

Held-out CHROMOSOME, strand+trinuc-matched negatives, base rate 0.0909:
  block             AUROC    top0.1%  top1%   top5%   top10%
  hairpin_only      0.5435   4.83x    2.545x  1.593x  1.333x
  sequence_only     0.6181   4.26x    2.750x  2.107x  1.827x
  hairpin+sequence  0.6205   4.96x    2.925x  2.143x  1.860x
  ALL random baselines 0.91-1.06 (flat)
  structure adds over sequence: +0.0024 AUROC
  sequence adds over structure: +0.0770 AUROC

INTERPRETATION (honest):
- Real BEYOND-MOTIF signal exists for endogenous A3A: 2.9x at top-1%, random 1.0x.
- Explicit hairpin features are largely REDUNDANT with a +/-10bp sequence model
  (+0.0024 AUROC). Hairpin-ness is itself a local-sequence property; the
  parameterisation is an interpretable summary, not extra information.
- At the extreme tail (top 0.1%) hairpin_only alone gets 4.83x vs combined 4.96x
  => the very top is essentially all structure.
- Magnitude is modest and the tail is small => this predicts HOTSPOTS, not burden.
  Consistent with the prior 0.28%-of-burden deterministic-core finding.
- This is the TRAINING-side validation only. The transfer test to editor clones
  is the actual question and awaits S2 alignment.

POSITIVE CONTROL ESTABLISHED: the harness detects real structural signal when it
exists, with flat random baselines. Any later null on the editor arms is now
interpretable rather than indistinguishable from a broken pipeline.

## RESULT — A3A vs A3B enzyme specificity (2026-08-10 00:18) — KEY RESULT
Third bug of the strand/convention family found LATENT and fixed: universe pos is
0-based, trainset pos is 1-based, same field name. Would have shifted every site
by 1 at the mpileup join (TCW -> CWN), silently. Added pos1 + CONVENTIONS.md.

Corrected strand-aware donor strata (background YTCA at TCW = 0.6057, DERIVED):
  A3A-like  ytca>=0.70 : 22 donors, 193,826 TCW muts, median 0.7208 = 1.190x bg
  A3B-like  ytca<=0.58 : 53 donors,  82,566 TCW muts, median 0.5083 = 0.839x bg
Both arms built at 80k positives + 10x strand+trinuc-matched negatives.
LEAK CHECKS PASS both arms (focal C 1.0000/1.0000; offsets -1,+1 identical).

HAIRPIN DOSE-RESPONSE, A3A-high vs A3B-like:
  stem>=5   A3A 1.234x   A3B 1.114x
  stem>=6   A3A 1.435x   A3B 1.203x
  stem>=7   A3A 1.709x   A3B 1.228x
  stem>=8   A3A 2.227x   A3B 1.221x     <- 1.82x separation
  pos_in_loop k=2: A3A 1.580x, A3B flat 1.123x
  all random baselines 0.85-1.25 (noisy only at n<150)

COMPOSITION CONFOUND EXCLUDED (strata differ: A3A 57.6% TCA, A3B 47.7% TCA):
  baseline hairpin rate in negatives is context-independent --
    TCA stem>=6 0.01672 / 0.01731 ; TCT stem>=6 0.01717 / 0.01695
  and A3A>A3B survives WITHIN each context:
    TCA stem>=8: A3A 2.721x (rand 1.01, n=280) vs A3B 1.300x (rand 1.02, n=116)
    TCT stem>=8: A3A 1.532x                    vs A3B 1.139x
    loop k=2 within TCA: A3A 1.698x vs A3B 1.160x
  => A3A effect concentrates in TCA; A3B has essentially NO hairpin dose-response.

MEANING: the enzyme-dichotomy premise (A3A-class has reference-predictable
structural preference that other deaminases lack) is VALIDATED on independent
cancer data. Previously it rested on a single AUROC 0.89-vs-0.51 comparison.
Gives the editor test a PRE-REGISTERED prediction: A3A-class editors should show
the steep TCA-hairpin dose-response; APOBEC1-class (YE1) and lamprey Lj-BE should
not; deaminase-free D10A/nCas9 should show nothing.
LIMITS: A3B stratum is not-A3A-like (low YTCA), not confirmed A3B. Effect is
hotspot-scale: 280 positives at stem>=8 of 80,000 = 0.35% of burden.

## CONTEXT-LENGTH SWEEP (in progress) — information-floor question
  +/-1bp AUROC=0.4989 (chance, as REQUIRED by trinuc matching = built-in null ctrl)
  +/-2bp 0.5862 (+0.0873)   <- the -2 base (YTCA/RTCA) carries almost everything
  +/-3bp 0.6021 (+0.0159)
  +/-5bp 0.6163 (+0.0142)
  +/-10bp 0.6181 (+0.0018)
  => saturating hard by +/-5bp. If this holds to +/-40, a DNA LM has little to add
     and A3A targeting is short-range sequence + something not in the reference.

## RESULT — A3A-vs-A3B enzyme claim: NULL (2026-08-10 08:40)
Pre-registered follow-up ran during the auth blackout. A3Awide = 102 donors
(ytca>=0.65) vs earlier 22 (ytca>=0.70). Leak checks pass.

Hairpin dose-response WEAKENS ~2x when the donor net widens:
  stem>=6  narrow(22d) 1.435x   wide(102d) 1.219x
  stem>=8  narrow(22d) 2.227x   wide(102d) 1.359x
  loop k=2 narrow      1.580x   wide       1.260x

DECISIVE donor-level test, n=97 (was 53):
  Spearman(enrichment, ytca)            = +0.079 p=0.44   NOT SIGNIFICANT
  Spearman(enrichment, burden)          = +0.087 p=0.39   NOT SIGNIFICANT
  PARTIAL enrichment~ytca  | burden     = +0.021 p=0.84
  PARTIAL enrichment~burden| ytca       = +0.021 p=0.83
  within burden tertiles rho(enr,ytca)  = +0.094 / -0.089 / +0.028 (all null)
  (ytca~burden collinearity +0.755)

VERDICT: hairpin enrichment does NOT track YTCA across 97 donors, in any burden
stratum. Both earlier correlations (+0.355/+0.365 at n=53) collapsed toward zero
as n grew => they were small-sample noise. The narrow-set 2.23x is best explained
as WINNER'S CURSE from selecting the extreme tail of a noisy statistic.
The A3A-vs-A3B enzyme-specificity claim is NULL. Two earlier characterisations in
this session (first 'validated', then 'suggestive') were both wrong.

SURVIVES UNCHANGED: the hairpin dose-response itself vs strand+trinuc-matched
negatives (flat random baselines, robust across all 10 GC deciles); the +/-5bp
context saturation; the positive control.

CONSEQUENCE FOR THE EDITOR TEST: it loses its confirmatory framing. The
pre-registered prediction (A3A-class shows TCA-hairpin dose-response, APOBEC1-class
does not) rested on this dichotomy validating. It did not. The editor test is now
EXPLORATORY. State that plainly; do not reframe after seeing the result.
NOTE: this does NOT refute the earlier A3A-hairpin result (AUROC 0.89 vs 0.51,
different data/analysis). What failed is reproducing the dichotomy via YTCA donor
stratification -- and YTCA/RTCA as an A3A-vs-A3B proxy is itself contested.

## MANDATORY ANALYSIS REQUIREMENT — germline filter (found 2026-08-10 09:1x)
QA of the calibrator's genome-wide counts (nCas9-clone1, 6 chroms, 77,154 alt>=2
sites) shows the raw call set is 35.5% GERMLINE-LIKE:
  VAF <0.05          0.074   noise/subclonal
  VAF 0.05-0.15      0.373   <- somatic/editing tier, where editing lives
  VAF 0.15-0.35      0.131
  VAF 0.35-0.65      0.132   germline het
  VAF >0.90          0.223   germline hom
  median VAF = 0.2148

WHY THIS WOULD HAVE BEEN FATAL, non-obviously: germline is SHARED between editor
and calibrator clones (same 293T parent), so it partly cancels in a ratio and
looks handled. But germline sites carry HIGH alt counts, so they dominate the top
of any alt-ranked list -- which is exactly our top-K% enrichment endpoint. The
editor signal (VAF 0.05-0.15) would be buried under germline homozygotes.

REQUIRED before ANY editor claim (mirrors the site-recurrence STEP-3 control that
salvaged the earlier Doman analysis):
  editor-specific site :=
      alt >= 2   in the editor clone
      AND alt == 0 AND cov >= 15 in Parent/background
      AND alt == 0 AND cov >= 8  in EVERY deaminase-free clone
        (nCas9-clone1/2 for PRJNA1042830; D10A-clone1/6/10 for PRJNA1006866)
This removes germline AND shared clonal APOBEC3 background in one step.
Do NOT rank by raw alt count without it.

Calibrator baseline (nCas9-clone1, genome-wide, 23 chroms):
  strand balance 1.0115 over 3.32M alt reads
  cov>=8 fraction 0.920-0.991 ; mean depth 17.6-32.2x
  alt>=2 rate 942-1994 /Mb, median 1140
Coverage confound at CHROMOSOME level: rho(rate, cov)=+0.258 p=0.23,
rho(rate, depth)=+0.312 p=0.15 -- both positive and n.s., i.e. no gross
detectability artifact. NOT a substitute for SITE-level coverage matching.

## CALIBRATOR BASELINE FOR HAIRPIN ENRICHMENT (2026-08-10 10:1x) — CRITICAL
s7_editor_test.py written and validated on a MUST-BE-NULL negative control:
nCas9-clone1 (deaminase-free) as pseudo-editor vs Parent. 23 chroms,
216,118,824 eligible background sites, 90,844 editor-specific sites.

  stem   enrich  RANDOM  n_spec
  >=4     0.938   0.995  13,848
  >=5     0.959   0.998   4,602
  >=6     1.051   0.971   1,601
  >=7     1.377   0.978     718
  >=8     1.625   0.994     296

A DEAMINASE-FREE clone shows 1.38-1.63x hairpin enrichment at high stem.
Random baselines flat (~0.98), so this is not a shuffling artifact.

TWO EXPLANATIONS, NOT YET SEPARATED:
 (a) endogenous APOBEC3A is active in these cells -- clone-private endogenous A3A
     mutations SHOULD show hairpin preference (same effect measured in PCAWG);
 (b) hairpin regions have worse mappability -> false calls in ANY clone.
Discriminating test (not yet run): compare VAF distribution and coverage of
hairpin vs non-hairpin specific sites. Artifact => low-VAF + coverage-dependent.

*** OPERATIONAL CONSEQUENCE ***
The deaminase-free baseline is ~1.4-1.6x, NOT 1.0. Every editor enrichment must
be quoted against THIS, never against unity. With the pre-registered >=3x
increment bar, an editor needs ~4.9x at stem>=8 to clear it.
Had the editor arm been run first, 1.6x could easily have been called an editor
effect. This surfaced ONLY because the analysis was validated on a case whose
answer was constrained in advance -- do that for every new analysis.
