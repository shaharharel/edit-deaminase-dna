# Pre-registered read-out for the editor test

Written 2026-08-22 ~04:50 UTC, **before D10A-clone6's counts exist** (58.7% aligned). The
whole point is that the decision rule is fixed while the numbers are still unknown. This
project has already had to retract one editor result that was framed after the fact.

The run fires unattended (`auto_advance2.sh`) and produces four things in order:
A the cross-study gate, B the deaminase-free validation, C the two editor arms unfiltered,
D the same two arms with the cross-clone recurrence filter on.

## Gate A — is Parent usable as the germline mask?

| concordance | conclusion | action |
|---|---|---|
| ≥ 0.97 | same lineage | proceed |
| 0.90–0.97 | caveat | proceed, but every downstream number carries the caveat in writing |
| < 0.90 | DIVERGED | **stop**; do not use Parent as mask; wait for P66-background (~14:10 UTC) |

Prior: Parent vs D10A-clone1 gave 0.9839 against a within-study reference of 0.9876, so ≥0.97
is expected. If it comes back below 0.90, nothing in C or D may be reported at all.

## Gate B — does the harness still return null on a null?

nCas9-clone2 vs nCas9-clone1, both deaminase-free, must be indistinguishable. Established
tonight with the banded estimator: MH **1.118/1.081** (stem 6), **1.383/1.380** (7),
**1.583/1.640** (8) — and verified against an independent re-implementation, 6/6 to three
decimals. If B does not reproduce those to within the ~12% noise floor, **the harness has
changed and C and D are not interpretable.** Report B first, always.

## Endpoint A — burden. The ≥3× bar lives here and only here.

Read **within matched coverage bands**, never pooled, and only in bands whose cov-ratio audit
sits inside 0.98–1.02. A band tagged `UNMATCHED DEPTH` does not carry the bar.

- **≥ 3× over the deaminase-free calibrator in depth-matched bands** → the first evidence
  tonight of an editor-attributable burden. Would still need clone10 (second calibrator,
  ~16:10 UTC) before leaving this document.
- **1.5–3×** → suggestive, not sufficient. Named as suggestive, not upgraded.
- **≈ 1× (0.88 is the clone floor)** → no detectable editor burden. This is what every arm has
  shown so far.

## Endpoint B — hairpin enrichment. **No 3× expectation, ever.**

Quoted against the deaminase-free calibrator's **1.28–1.63×**, never against 1.0. An
A3A-family editor that targets like endogenous A3A should MATCH the calibrator, not exceed it.

- **Editor MH clearly above the calibrator's** at the same stem, outside the noise floor, in
  both arms → a real difference worth pursuing.
- **Editor at or below the calibrator** → consistent with everything measured tonight:
  substrate, not enzyme. Lj-BE (lamprey CDA1, a different deaminase family) sits at the
  calibrator too, and that is the strongest argument that this ~1.1–1.6× is a property of
  which sites get called in clonal WGS.
- **n_hp below ~100 in a band** → report the number with its n and draw nothing from it. The
  A3A-vs-A3B claim died at exactly this size.

## C vs D — what the recurrence filter tells us

Two independent A3A-Y130F clones share 3,000 editor-specific sites against 1.8 expected.
Splitting on that sharing gave 1.326× shared vs 1.026× private at stem 8.

- **C ≈ D** → the arm's signal is clone-private and survives the control.
- **C > D materially** → the signal lives in the shared sites, i.e. the artefact class, and
  the arm is downgraded exactly as it was before. Note that eA3A behaves the *opposite* way
  (private 1.371× vs all-three 1.110×), so this is a real discriminator, not a formality.

## What cannot be concluded regardless of the numbers

- Nothing about **A3A vs A3B**. That dichotomy failed to validate and is not being retested.
- Nothing **confirmatory**. The editor test is EXPLORATORY because the dichotomy failed; a
  positive result is a hypothesis for a designed experiment, not a finding.
- **haA3A (Y130G/VA) was engineered for near-background off-target — expect NO enrichment.**
  If it shows one, suspect the pipeline before the biology.
- Nothing from a single clone. clone7 (~15:50) and clone10 (~16:10) are what turn any of this
  from n=2 into n=3 with a second calibrator.
