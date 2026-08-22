# READ THIS FIRST — state as of 2026-08-22 03:00 UTC

STATUS.md is now ~3,700 lines of append-only tick log. This page is the current state; the
log is the evidence trail behind it. Written to be read in two minutes.

## What finished

**The DeaminaFormer architecture search is COMPLETE and decisive.** A pretrained DNA
foundation model adds nothing at the operating point.

| input block | MLP head | GB head |
|---|---|---|
| hairpin geometry (6 feat) | 4.881× | 4.940× |
| hairpin + thermodynamic (16 feat) | 4.655× | **5.036×** |
| NT-v2 embeddings alone (1024 dim) | 1.250× | **1.048×** |
| NT-v2 + hairpin geometry | 1.488× | 4.857× |
| everything fused | 2.333× | 4.774× |
| GB one-hot + hairpin (prior baseline) | — | 5.280× |

random baseline 0.917× · arithmetic ceiling 11.000× · n=924 · 5 folds by held-out chromosome

The embeddings are EMPTY, not diluted — gradient boosting can ignore noise features and still
gets 1.048×. The MLP's collapse when embeddings were added was the HEAD (GB: 4.940 → 4.857).
Nothing beats hairpin geometry, including a Conv1d motif scanner over ±200 bp.

**Recommended model: 16 features, gradient-boosted.** 5.036× at top 0.1%, ECE 0.00125, and
calibrated where it operates (top 0.1%: predicts 0.4822, observes 0.4782). Auditable and
cheap — a better outcome for a regulatory gate than a foundation model.

**Published**: HANDOFF §14 on branch `phase7-a3a-hairpin` (commit ec2e800), 55 scripts and
three result JSONs. 13/13 table values verified against the artifacts.

## What is still running

| sample | node | done | ETA (UTC) | unblocks |
|---|---|---|---|---|
| P66-D10A-clone6 | A | 53.0% | ~10:10 | **the editor test** — clean calibrator |
| P66-Lj-BE-clone5 | B | 83.0% | ~04:25 | cross-family control |
| P66-Lj-BE-clone12 | B | 58.2% | ~07:50 | cross-family control |
| P66-background | A | 42.0% | ~14:10 | matched germline mask |
| P66-A3A-Y130F-clone7 | A | 38.6% | ~15:50 | stem-8 n beyond 96 |
| P66-D10A-clone10 | A | 38.1% | ~16:10 | second clean calibrator |

Both nodes have an **auto-advance driver armed**: the analyses fire the moment counts complete,
without waiting for a monitor tick. Results land in `logs/auto_advance.log` (node A) and
`logs/auto_advance_B.log` (node B). `ai-gpu` is STOPPED; its scratch disk persists.

## Where the editor question stands — NO CLAIM, and the framing has not moved

Every arm sits at or below the deaminase-free calibrator (1.28–1.63× at stem 8). Tonight:

- **Banded endpoint B is now the method of record.** `alt≥2` is a depth-dependent VAF
  threshold, so arms of unequal depth were never comparable. A VAF floor FAILED its own
  validation (two deaminase-free clones diverged 1.438 vs 1.159). Matched coverage bands
  passed it (1.118/1.081, 1.383/1.380, 1.583/1.640) and were verified against an independent
  re-implementation, 6/6 values to three decimals.
- **RETRACTED**: the stem-6 between-arm difference (A3A-Y130F 0.877 vs eA3A 1.061,
  p=2.6×10⁻⁶). Controlling VAF reverses the direction. It measured VAF composition.
- **Cross-family, depth-matched**: Lj-BE (lamprey CDA1) 1.138/1.308/1.569 vs eA3A
  1.129/1.210/1.324 — indistinguishable, both at the calibrator. Under the pre-registration
  that is the "substrate, not enzyme" branch. Still n=1 Lj-BE clone.
- **The eA3A enrichment is clone-PRIVATE** (1.371× private vs 1.110× all-three at stem 8) —
  the opposite of A3A-Y130F, whose signal was carried entirely by shared sites and was
  downgraded for it.

## Three things that will bite whoever picks this up

1. **The empirical noise floor is ~12%** at top 0.1%. The random baseline, which should read
   1.000×, is 0.973 ± 0.062 across 18 blocks. Differences under 12% are not resolvable — this
   covers 5.036 vs 5.280, and it already retired one claim tonight.
2. **Calibration is to the sampled 1-in-11 base rate**, not genome-wide prevalence. Re-base
   any deployed threshold before quoting a risk number.
3. **Tier A/B recall — the spec's headline metric — is blocked on DATA**, not analysis:
   COSMIC Tier 1, ClinGen HI Level 3, DepMap common essentials. None on either node.

## Instruments you should trust, and one you should not

- `progress.sh` (both nodes) — percentages from live fd offsets ÷ real file size. Authoritative.
- Its `CUR_READS` column **double-counts any sample re-aligned with the same thread count**
  (clone6 reads 2748M implied against a true ~991M). The `IMPLIED` column exists to expose
  exactly that; compare it across samples with similar file sizes.
- The cron's read-count snippet uses assumed targets and cannot see restarts. Use `progress.sh`.
