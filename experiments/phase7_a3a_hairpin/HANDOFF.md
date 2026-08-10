# HANDOFF — Phase 7: A3A hairpin / PCAWG → base-editor transfer

_Written 2026-08-10 ~10:40 UTC. Read this first when resuming in a new session._

## 0. THE ONE THING THAT BREAKS ON SESSION SWITCH

**Node-side work keeps running. Claude-side monitoring does not.**

- All compute on `ai-chem` runs as **systemd units** — aligners, pileup driver,
  reaper. These survive any Claude session ending. Nothing stops.
- The **monitor and QA crons are session-only** (in-memory, die with the session).
  A new session MUST recreate them or the run proceeds unwatched.

Recreate both immediately on resume — the exact prompts are in §5.

## 1. Where the work lives

| what | where |
|---|---|
| all scripts | `experiments/phase7_a3a_hairpin/scripts/` (this repo) AND `/mnt/data/a3a/` on ai-chem |
| running status, all results, bug post-mortems | `/mnt/data/a3a/STATUS.md` (copy in scripts/) |
| coordinate conventions — READ BEFORE JOINING ARRAYS | `/mnt/data/a3a/CONVENTIONS.md` |
| counts, universes, trainsets | `/mnt/data/a3a/feat/` (node only, ~100 GB) |
| logs | `/mnt/data/a3a/logs/` |
| idempotency flags | `/mnt/data/a3a/flags/` |

Node: `gcloud compute ssh ai-chem --zone us-east1-b --tunnel-through-iap --command="..."`
Python: `~/miniconda3/envs/apobec/bin/python`

## 2. What this phase is testing

Train a sequence/structure classifier on **endogenous APOBEC3A mutations from
cancer genomes (PCAWG)**, then ask whether it transfers to **base-editor
off-targets**. The governing idea: hairpins make ssDNA exposure
**sequence-determined**, which is the one mechanism that can transfer. Every
cell-state ssDNA proxy tried in earlier phases was null.

## 3. Results so far

**Stands:**
- **Hairpin dose-response is real.** Enrichment 1.05× (stem≥3) → 2.32× (stem≥9)
  for APOBEC mutations vs strand+trinucleotide-matched negatives, flat random
  baselines, survives all ten GC deciles and complexity matching.
- **Sequence context saturates at ±5 bp.** AUROC 0.4989 (±1) → 0.6163 (±5) →
  0.6185 (±15) → 0.6178 (±40). Thirty extra bp buy +0.0007. **A DNA LM has
  nothing to find** — this is the information floor, measured directly.
- **Reproduced independently in HEK293T clones.** The deaminase-free calibrator
  shows 1.38–1.63× hairpin enrichment at stem≥7/≥8, and a VAF/coverage test
  showed it is genuine endogenous A3A, not mappability artifact.

**Dead:**
- **A3A-vs-A3B enzyme dichotomy is NULL.** ρ(enrichment, YTCA) = +0.355 at n=53
  collapsed to **+0.079 (p=0.44) at n=97**; partials +0.021/+0.021; null in every
  burden tertile. The earlier 22-donor 2.23× was **winner's curse** from selecting
  the extreme tail of a noisy statistic. I called this "validated" and then
  "suggestive" before it resolved to null — both were wrong.

**Consequence:** the editor test is **EXPLORATORY, not confirmatory**. Do not
restore confirmatory framing after seeing results.

## 4. Non-negotiable discipline (all learned the hard way here)

1. **Negatives matched on trinucleotide context AND strand.** Unmatched → a
   meaningless AUROC that is pure motif-learning.
2. **Random baseline beside every enrichment.** A number equal to `1/base_rate`
   is the arithmetic ceiling = leakage, not skill (we hit exactly 11.000×).
3. **Derive backgrounds, never assume.** YTCA background at TCW is **0.6057**,
   not 0.5. Assuming 0.5 made a 0.55 threshold silently select noise.
4. **Report n with every effect** and check whether it strengthens or weakens as
   n grows. Small-n illusions killed the A3A/A3B claim.
5. **Enrichment is the endpoint, not AUROC.** They dissociate; trust enrichment.
6. **Germline filter is mandatory** — 35.5% of raw alt≥2 calls are germline-like,
   and germline carries high alt counts so it dominates top-K rankings.
   Editor-specific := `alt≥2 & cov≥8` in editor, `alt==0 & cov≥15` in
   Parent/background, `alt==0 & cov≥8` in EVERY deaminase-free control.
7. **Calibrator baseline is ~1.4–1.6×, NOT 1.0.** Pre-registered bar is a **≥3×
   increment over deaminase-free**, so an editor needs ~4.9× at stem≥8.
8. **Validate every new analysis on a case whose answer is constrained in
   advance.** The calibrator baseline was found only this way.
9. **Subagents do not work in this project's sessions** — all six spawned went
   idle without returning content. Do QA inline.

## 5. Recreate the crons on resume

**Monitor, every 10 min** (`4,14,24,34,44,54 * * * *`) and **QA, every 30 min**
(`9,39 * * * *`). Full prompt text is preserved in `cron_prompts.md` next to this
file. The QA prompt must say *do the QA yourself, do not spawn subagents*.

## 6. Live state at handoff

Aligning on ai-chem (measured rates, not estimates):

| sample | progress | ETA |
|---|---|---|
| nCas9-clone2 (calibrator) | 98% | minutes |
| **P66-A3A-Y130F-clone2 (EDITOR)** | 93.0M / 825.9M reads | ~7.6 h |
| **P66-D10A-clone1 (matched calibrator)** | 40.2M / 749.0M | ~9.7 h |
| P66-background (4 threads, laggard) | 240.3M / 942.2M | ~31 h |

Complete: Parent + nCas9-clone1 counts (23/23 chromosomes each).
**First editor result ≈ 10 h from handoff.** Full 19-sample queue ≈ 2.5 days.

## 7. Open decisions (I recommended, user had not answered)

1. **Restart P66-background at 9 threads** — costs 10.5 h spent, finishes in ~16 h
   vs 31 h remaining. Net saving ~14 h. Capacity frees as nCas9-clone2 completes.
2. **Redirect w3 to the queue** when it finishes nCas9-clone2, otherwise it starts
   `YE1-clone1` (hardcoded, expected-null). w1 and w2 were already redirected this
   way; `queue.txt` is ordered by information value and drops nothing.
3. **Cross-study germline masking is untested.** Parent (PRJNA1042830) as germline
   mask for P66 (PRJNA1006866) samples assumes a shared HEK293T lineage. Our
   98.76% homozygous concordance was measured *within* PRJNA1042830. Verify
   cross-study concordance as soon as a P66 sample has counts.

## 8. Known bugs found and fixed (all the same family)

A convention correct locally, then consumed downstream as if universal. None
crashed; each would have returned a confident wrong number.

1. Negatives sampled plus-strand-only → exactly-ceiling 11.000× enrichment.
2. `-2` base via a `ref=='C'` branch → YTCA/RTCA and the donor split corrupted.
3. universe `pos` 0-based vs trainset `pos` 1-based, same field name (caught latent).
4. `s6_pileup` counting only uppercase `T` → would drop half the reads
   strand-asymmetrically (the old Selict 21.8:1 mechanism).
5. `pos1` missing from 22 of 23 universes — I patched the artifact for bug 3 but
   not the generator, so it recurred 22×. **Fix generators, not artifacts.**
6. The reference-base guard was **vacuous** — `mpileup -f REF` reports the base
   from the FASTA, not the reads, so it compared a FASTA to itself. Replaced with
   a BAM `@SQ`-vs-`.fai` build check. **A guard that has never fired is an
   assumption with a comment attached.**
