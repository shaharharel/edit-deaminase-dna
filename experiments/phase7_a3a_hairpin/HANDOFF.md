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

**Stands (PCAWG / training side):**
- **Hairpin dose-response is real.** 1.05x (stem>=3) -> 2.32x (stem>=9) for APOBEC
  mutations vs strand+trinucleotide-matched negatives, flat random baselines,
  survives all ten GC deciles and complexity matching.
- **Sequence context saturates at +/-5 bp.** AUROC 0.4989 (+/-1) -> 0.6163 (+/-5)
  -> 0.6185 (+/-15) -> 0.6178 (+/-40). **A DNA LM has nothing to find** -- this is
  the information floor, measured rather than inferred.

**Stands (HEK293T clones -- independent reproduction):**
- The deaminase-free calibrator shows **1.28-1.63x** hairpin enrichment at
  stem>=7/>=8, **replicated across two independent clones** (1.377/1.625 vs
  1.280/1.531, agreeing to ~6%).
- It is **endogenous A3A biology, not artifact**: hairpin-specific sites match
  non-hairpin on VAF (0.0800 vs 0.0769) and coverage (27.0 vs 28.0), and the
  effect survives coverage stratification (1.31-1.57x across 8-15x/15-25x/
  25-35x/35x+ bins, n=3.6k-36k per bin).
- **Statistically established, not just observed**: with a 2000-draw null
  distribution (s7b_stats.py), observed CIs sit entirely above null CIs at
  stem>=6/7/8 (p<0.0001 at >=7 and >=8).
- **Survives within trinucleotide context** -- TCA 1.144/1.345/1.660 and
  TCT 1.042/1.207/1.379 at stem>=6/7/8. Not composition. Effect is STRONGER in
  TCA, the same asymmetry seen in PCAWG.
- **98.76% germline-homozygous concordance** between two independently processed
  genomes, while the low-VAF (editing) tier is only 6.8% shared -- i.e. clone-
  private, which is the precondition for detecting editor-specific events.

**Dead:**
- **A3A-vs-A3B enzyme dichotomy is NULL.** rho(enrichment, YTCA) = +0.355 at n=53
  collapsed to **+0.079 (p=0.44) at n=97**; partials +0.021/+0.021; null in every
  burden tertile. The 22-donor 2.23x was **winner's curse**. I called this
  "validated", then "suggestive", before it resolved to null -- both were wrong.

**Consequence:** the editor test is **EXPLORATORY, not confirmatory**.

## 3b. THE BAR, and the traps around it

**Editor must reach ~4.6-4.9x at stem>=8** -- i.e. >=3x over the replicated
deaminase-free baseline of 1.28-1.63x. **Never quote an editor number against 1.0.**

**TRAP 1 -- configuration dependence (3x effect, dangerous).** A BULK sample
(Parent) scored against clonal controls gives **4.941x** at stem>=8 (CI 4.66-5.22)
-- essentially ON the bar, from population structure alone. Bulk pools endogenous
A3A across many lineages so recurrent hairpin hotspots accumulate; one clone
carries one lineage. NEVER compare enrichments across sample configurations; the
calibrator must match the editor in clonality.

**TRAP 2 -- control count (3% effect, ignorable).** 1 vs 2 controls gives
1.574 vs 1.531 at stem>=8, well inside overlapping CIs. So using more controls in
the editor test than the calibrator had does NOT require re-deriving the baseline.

**Known conservative bias:** the control mask preferentially removes hairpins
(retention 0.9803 at stem>=8, 0.9737 at stem>=9, vs 0.9896 overall). 1-2% against
a >=3x requirement -- it can only understate a real effect, never invent one.

**Reporting requirement:** stratify by trinucleotide context, not just pooled --
specific sites run 52.2% TCA vs 47.0% background, so a mix shift could imitate an
editor effect.

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

## 6. Live state (updated 2026-08-10 ~15:50 UTC)

Node uptime ~16.5 h, never preempted. ~514 GB free. Seven systemd units active.

| sample | progress | role |
|---|---|---|
| **P66-A3A-Y130F-clone2** | 534.6M / 825.9M (65%) | **EDITOR** |
| **P66-D10A-clone1** | 414.0M / 749.0M (55%) | **matched calibrator** |
| P66-A3A-Y130F-clone5 | 286.8M / 833.1M | editor replicate |
| P66-background (4 threads) | 393.6M / 942.2M | matched parent, laggard |

**Complete with 23/23 chromosome counts:** Parent, nCas9-clone1, nCas9-clone2
(all three calibrators). Their BAMs were freed after pileup; counts retained.

Running units: `a3a-s2-q1..q4` (aligners), `a3a-s6-driver2` (pileups + selective
BAM retention), `a3a-reaper` (stale-claim release), `a3a-q1-watchdog`.

**First editor result** needs A3A-Y130F-clone2 + D10A-clone1 both piled up.
Per-sample rates under 4-way contention are ~60-70M reads/h, so roughly 5-7 h
from this timestamp, plus ~30 min of pileup each. Full 19-sample queue ~2 days.

## 7. Decisions taken (previously open)

1. **w1, w2 and w3 were all redirected** from hardcoded expected-null arms
   (Y130G / VA / YE1) onto the queue, each stopped during download so nothing
   aligned was discarded. All four alignment slots now serve the decisive
   contrast. `queue.txt` is ordered by information value; nothing was dropped.
2. **P66-background was NOT restarted at 9 threads.** I recommended this earlier
   on the reasoning that it would save ~13 h -- that reasoning was WRONG. Threads
   are zero-sum on a 32-core box, so reallocating would have taken them from other
   samples without changing total work. The real lever was queue order, not
   thread count.
3. **`a3a-q1-watchdog` armed.** q1 still runs the OLD queue script (FTP URLs, no
   FAILED guard). ENA's FTP now 404s for PRJNA1006866 paths while HTTPS works, and
   the old script re-claims instantly on failure -- it spun 9x in 11 s when this
   first happened. The watchdog retires q1 the moment P66-background completes (or
   if BLOCKED.md exceeds 20 lines) and starts a replacement on `s2_queue_v2.sh`.
   **This fires after the session is likely to have ended -- do not remove it.**
4. **Pileup driver v2 retains editor-arm BAMs** (`A3A-Y130F`) after their 23
   chromosomes, freeing controls as before, guarded by a 150 GB floor. Counts
   cover most follow-ups but not read-level inspection, and regenerating a BAM
   is ~8 h against 21 GB of disk.

## 7b. Still untested

- **Cross-study germline concordance.** The first editor test will likely use
  Parent (PRJNA1042830) as germline mask for P66 (PRJNA1006866) samples, because
  P66-background lags. Our 98.76% homozygous concordance was measured WITHIN
  PRJNA1042830. 293T sublines drift between labs. **Verify as soon as any P66
  sample has counts** -- same check, Parent vs a P66 control.

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
