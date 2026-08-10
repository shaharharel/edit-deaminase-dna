# HANDOFF — Phase 7: A3A hairpin / PCAWG → base-editor transfer

_Last updated 2026-08-10 ~19:00 UTC. Read §0 and §1 before doing anything._

---

## 0. THE ONE THING THAT BREAKS ON SESSION SWITCH

**Node-side work keeps running. Claude-side monitoring does not.**

- Everything on `ai-chem` runs as **systemd units** — aligners, pileup driver,
  reaper, watchdog. These survive any session ending. Nothing stops.
- The **monitor and QA crons are session-only** (in-memory). They die with the
  session and MUST be recreated, or the run proceeds unwatched.

Verbatim prompts: `cron_prompts.md` (monitor every 10 min, QA every 30 min).
The QA prompt must keep the line *"do the QA yourself, do NOT spawn subagents"* —
every subagent spawned in this project went idle without returning content.

**Also armed and load-bearing:** `a3a-q1-watchdog` fires **after this session is
likely to have ended**. Do not remove it (see §9.3).

---

## 1. WHY THIS PHASE EXISTS (context a new session will not have)

Earlier phases tried to predict **guide-independent** base-editor off-targets from
the reference genome and failed everywhere: Doman BE4/YE1 WGS, Lei Detect-seq,
McGrath AncBE4max, Selict ABE8e — all null under model ranking. The decomposition
showed why: apparent enrichment was **coverage artifact + endogenous clonal
APOBEC3 background**, and the localising variable is **τ_ssDNA** — whether a base
happened to be single-stranded at that moment — which is *not in the reference*.

The advisor's framing ("this is RNA") was that the whole approach was the RNA
editing-index paradigm transplanted onto DNA, where it does not hold: RNA editing
is steady-state and re-drawn continuously, DNA editing is one-shot and lineage-
recorded, and the DNA control is not flat because endogenous APOBEC3 targets the
same ssDNA.

**Phase 7 is the one mechanism that escapes that floor.** A hairpin makes ssDNA
exposure **sequence-determined** rather than cell-state stochastic. If A3A targets
hairpins, that preference IS computable from the reference and should transfer.
Train on endogenous A3A mutations in cancer genomes (PCAWG), then test transfer to
base-editor off-targets in clonal WGS.

---

## 2. PIPELINE ARCHITECTURE

Stages are independent and idempotent via flag files in `/mnt/data/a3a/flags/`.
Node: `gcloud compute ssh ai-chem --zone us-east1-b --tunnel-through-iap --command="..."`
Python: `~/miniconda3/envs/apobec/bin/python`. Aligner env: `~/miniconda3/envs/bio/bin`.

```
S1  PCAWG training set        S4  scoring universe        S2  clone WGS
    pcawg_build.py                s4_universe.py              s2_queue_v2.sh
    -> parse 21.6M SNVs           -> every TCW C, both        -> HTTPS download,
    fix_tetra.py                     strands, per chrom          gzip -t verify,
    -> strand-aware -2 base       -> hairpin feats            -> bwa mem, sort,
       donor_ranking_v2.tsv          precomputed                 index
    pcawg_stage2b.py              (23 chroms DONE)           reaper.sh, q1_watchdog.sh
    -> 100k pos + 1M matched neg
    build_A3A/A3B/A3Awide.py
              |                                 |                    |
              v                                 v                    v
S3  classifier                            S6  pileup: s6_pileup.py per (sample,chrom)
    s3_train_v2.py                            orchestrated by s6_driver_v2.sh
    -> hairpin/seq ablation                   -> counts_<sample>_chr<N>.npz
S5  s5_ctxsweep.py                               (cov, alt, alt_fwd, alt_rev)
    -> context saturation                     -> frees control BAMs, RETAINS
                                                 A3A-Y130F BAMs
                                                          |
                                                          v
                          S7  EDITOR TEST -- use s7c_editor.py (NOT s7/s7b)
                              -> endpoint A burden (coverage-matched)
                              -> endpoint B hairpin enrichment (null dist, by context)
```

**Script status — use the right version:**

| use | do NOT use | why |
|---|---|---|
| `s7c_editor.py` | `s7_editor_test.py`, `s7b_stats.py` | only s7c does coverage matching; s7 also lacks the null distribution |
| `s2_queue_v2.sh` | `s2_queue.sh`, `s2_wgs.sh` | v2 uses HTTPS (ENA FTP 404s on PRJNA1006866) and has a FAILED-flag backoff; v1 spin-loops |
| `s6_driver_v2.sh` | `s6_driver.sh` | v2 retains editor BAMs |
| `pcawg_stage2b.py` | `pcawg_stage2.py` | stage2 had the plus-strand-only negatives bug |
| `a3a_donor_ranking_v2.tsv` | `a3a_donor_ranking.tsv` | v1 YTCA values corrupted by the −2 base bug |

---

## 3. DATA INVENTORY (`ai-chem:/mnt/data/a3a/`)

| path | contents |
|---|---|
| `pcawg/pcawg_icgc_public.maf.gz` | 925 MB, PCAWG open tier, **no DACO needed** (`s3://icgc25k-open`, endpoint `https://object.genomeinformatics.org`) |
| `pcawg/snvs.npz` | 113 MB, 21,606,842 parsed SNVs, pyrimidine-oriented, hg19 |
| `pcawg/minus2.npy` | 86 MB, −2 base — **v1, CORRUPTED**; use `fix_tetra.py` output |
| `pcawg/a3a_donor_ranking_v2.tsv` | strand-aware per-donor tcw_frac / n_tcw / ytca_frac |
| `feat/universe_chr*.npz` (23) | every TCW C both strands + hairpin features. Fields: `pos` (0-based), `pos1` (1-based, **join on this**), `strand`, `tri` (0=TCA,1=TCT), `stem`, `loop`, `kpos`, `gc_pairs`, `hp_score`, `local_gc` |
| `feat/universe_chr*.bed` (23) | 5.3 GB total, for `mpileup -l` |
| `feat/a3a_trainset_v2.npz` | 100k pos + 1M strand+trinuc-matched neg, incl. `win` (81 bp strand-oriented windows) |
| `feat/a3a_trainset_A3A/A3B/A3Awide.npz` | donor strata |
| `feat/counts_<sample>_chr<N>.npz` | per-site `cov`, `alt`, `alt_fwd`, `alt_rev`; ~13 MB each. **3 samples × 23 chroms done** |
| `CONVENTIONS.md` | **read before joining any two arrays** |

Reference `/mnt/data/ref/hg19.fa` (+ bwa index). **Everything is hg19/GRCh37**, matching PCAWG.

---

## 4. RESULTS

### Stands — PCAWG training side
- **Hairpin dose-response real**: 1.05× (stem≥3) → 2.32× (stem≥9) vs
  strand+trinucleotide-matched negatives, flat random baselines, survives all ten
  GC deciles and complexity matching.
- **Sequence context saturates at ±5 bp**: AUROC 0.4989 (±1) → 0.6163 (±5) →
  0.6185 (±15) → 0.6178 (±40). **A DNA LM has nothing to find.** Information floor
  measured directly rather than inferred from failed models.
- Structure adds only **+0.0024 AUROC** over a ±10 bp sequence model — hairpin-ness
  is itself a local-sequence property. But the extreme tail is *all* structure
  (hairpin-only 4.83× vs combined 4.96× at top 0.1%).

### Stands — HEK293T clones, independent reproduction
- Deaminase-free calibrator shows **1.28–1.63×** at stem≥7/≥8, **replicated across
  two independent clones** (1.377/1.625 vs 1.280/1.531, ~6% agreement).
- **Endogenous A3A, not artifact**: hairpin-specific sites match non-hairpin on VAF
  (0.0800 vs 0.0769) and coverage (27.0 vs 28.0); survives coverage stratification
  (1.31–1.57× across 8–15/15–25/25–35/35×+ bands, n=3.6k–36k).
- **Statistically established**: 2000-draw null; observed CIs entirely above null
  CIs at stem≥6/7/8 (p<0.0001 at ≥7 and ≥8).
- **Survives within trinucleotide context**: TCA 1.144/1.345/1.660, TCT
  1.042/1.207/1.379 at stem≥6/7/8. Stronger in TCA — same asymmetry as PCAWG.
- **98.76% germline-hom concordance** between independently processed genomes,
  while the low-VAF editing tier is only **6.8%** shared (clone-private, the
  precondition for detecting editor-specific events).

### Dead
- **A3A-vs-A3B enzyme dichotomy is NULL.** ρ(enrichment, YTCA) +0.355 at n=53 →
  **+0.079 (p=0.44) at n=97**; partials +0.021/+0.021; null in every burden
  tertile. The 22-donor 2.23× was **winner's curse**. Called "validated" then
  "suggestive" before it resolved to null — both were wrong.
- **Consequence: the editor test is EXPLORATORY, not confirmatory.** Do not
  restore confirmatory framing after seeing results.

---

## 5. THE BAR, AND THE TRAPS AROUND IT

**Two endpoints, different confounds. Report both.**

**(A) BURDEN** — editor-specific site count. **The ≥3× bar belongs here.**
Calibrator reference 90,844 sites (1 control) / 73,512 (2 controls), so a passing
editor shows ~220–270k. **Must be coverage-matched**: P(alt≥2 | VAF 0.08) is 0.483
at 20× and 0.841 at 40×, so depth alone inflates burden **1.74×**, and the P66
editor samples are ~1.4× deeper than the calibrators.

**(B) HAIRPIN ENRICHMENT** — a within-sample ratio, largely depth-immune. Quote
against the calibrator's **1.28–1.63×**, *not* against 1.0, and with **no 3×
expectation**: an A3A-family editor that targets like endogenous A3A should MATCH
the calibrator. ≈1.0 refutes the hypothesis; ≫ calibrator would be a strong new
claim.

**Informative pattern = A high AND B near the calibrator.** B large with A flat is
more likely artifact than editing.

**Noise floor (measured):** two deaminase-free clones differ **~12% in burden,
~5% in enrichment**. The 3× bar sits far above that.

**TRAP 1 — configuration dependence (3×, dangerous).** A **bulk** sample scored
against clonal controls gives **4.941×** at stem≥8 (CI 4.66–5.22) — essentially ON
the bar, from population structure alone. Bulk pools endogenous A3A across many
lineages so recurrent hairpin hotspots accumulate. **Never compare enrichments
across sample configurations**; the calibrator must match the editor in clonality.

**TRAP 2 — control count (3%, ignorable).** 1 vs 2 controls: 1.574 vs 1.531. Using
more controls than the calibrator had does not require re-deriving the baseline.

**Known conservative bias:** the control mask preferentially removes hairpins
(retention 0.9803 at stem≥8, 0.9737 at stem≥9, vs 0.9896 overall) — 1–2% against a
≥3× requirement. Can only understate a real effect.

**Reporting requirement:** stratify by trinucleotide context. Specific sites run
52.2% TCA vs 47.0% background; a mix shift could imitate an editor effect.

---

## 6. NON-NEGOTIABLE DISCIPLINE

1. **Negatives matched on trinucleotide context AND strand.** Unmatched → a
   meaningless AUROC that is pure motif-learning.
2. **Random baseline beside every enrichment.** A value equal to `1/base_rate` is
   the arithmetic ceiling = leakage, not skill (we hit exactly 11.000×).
3. **Derive backgrounds, never assume.** YTCA background at TCW is **0.6057**, not
   0.5. Assuming 0.5 made a 0.55 threshold silently select noise.
4. **Report n with every effect**; check whether it strengthens or weakens as n
   grows. Small-n illusions killed the A3A/A3B claim.
5. **Enrichment is the endpoint, not AUROC.** They dissociate.
6. **Germline filter is mandatory** — 35.5% of raw alt≥2 calls are germline-like
   and carry high alt counts, so they dominate top-K rankings.
   `alt≥2 & cov≥8` in editor; `alt==0 & cov≥15` in Parent/background;
   `alt==0 & cov≥8` in EVERY deaminase-free control.
7. **Calibrator baseline is ~1.4×, NOT 1.0.**
8. **Validate every new analysis on a case whose answer is constrained in
   advance.** The calibrator baseline and the noise floor were both found this way.
9. **Subagents do not work here.** Do QA inline.
10. **Fix generators, not artifacts.** Bug 3 was "fixed" in one file and recurred
    22 times because the generator was untouched.
11. **A guard that has never fired is an assumption with a comment attached.**
    Bug 6 was a vacuous check that compared a FASTA to itself.

---

## 7. HOW TO RUN THE EDITOR TEST

Once an editor and its matched calibrator both have 23/23 counts:

```bash
~/miniconda3/envs/apobec/bin/python /mnt/data/a3a/s7c_editor.py \
    P66-A3A-Y130F-clone2 \
    P66-D10A-clone1 \
    P66-D10A-clone6,P66-D10A-clone10 \
    P66-background
# args: <editor> <calibrator> <extra_controls_csv|''> <parent>
```

**Before trusting any output, re-run in validation mode** — two deaminase-free
clones as editor/calibrator. Burden ratio must be ~1.0 per band (measured 0.88,
the clone-to-clone floor) and enrichment must reproduce ~1.5×:

```bash
... s7c_editor.py nCas9-clone2 nCas9-clone1 '' Parent
```

If `P66-background` is not yet aligned, `Parent` (PRJNA1042830) can serve as the
germline mask — but see §11: cross-study concordance is **UNTESTED**.

---

## 8. LIVE STATE (2026-08-10 ~19:00 UTC)

Node uptime ~20 h, never preempted. ~495 GB free. Seven units active.

**Complete, 23/23 counts:** `Parent`, `nCas9-clone1`, `nCas9-clone2` (all three
calibrators). BAMs freed; counts retained.

| sample | progress | role |
|---|---|---|
| P66-A3A-Y130F-clone2 | 723.6M / 825.9M (88%) | **EDITOR** |
| P66-D10A-clone1 | 627.6M / 749.0M (84%) | **matched calibrator** |
| P66-A3A-Y130F-clone5 | 556.8M / 833.1M | editor replicate |
| P66-background | 482.7M / 942.2M (4 threads) | matched parent, laggard |

**Queue: 19 samples, none of the P66 set finished yet.** Order by information
value: A3A-Y130F ×3 interleaved with D10A ×3, then eA3A-RL1 ×3, Lj-BE ×3, then the
expected-null PRJNA1042830 arms (Y130G, VA, YE1), P66-background last.

**First editor result ≈ 1–2 h** from this timestamp. Full queue ≈ 2 days.

---

## 9. RUNNING UNITS — what they are and why

1. **`a3a-s2-q1..q4`** — queue-driven aligners. q1 runs the OLD script (see 9.3).
2. **`a3a-s6-driver2`** — sweeps any sample with `S2_*_DONE` through 23 chromosome
   pileups at 4-way parallelism, then frees the BAM — **except** samples matching
   `A3A-Y130F`, retained (150 GB floor) because read-level follow-up on the editor
   arm would otherwise cost ~8 h to regenerate.
3. **`a3a-q1-watchdog`** — **do not remove.** q1 runs the old queue script with FTP
   URLs and no FAILED guard. ENA's FTP 404s on PRJNA1006866 paths (HTTPS works),
   and the old script re-claims instantly on failure — it spun 9× in 11 s when this
   first happened. The watchdog retires q1 when P66-background completes (or if
   `BLOCKED.md` exceeds 20 lines) and starts a replacement on `s2_queue_v2.sh`.
   **It fires after this session is likely to have ended.**
4. **`a3a-reaper`** — releases `S2_*_CLAIMED` flags older than 4 h **only when no
   live process matches the sample**. Age alone is insufficient: a healthy 6.5 h
   job would otherwise be released and duplicated by a second worker.

**Diagnosing a "stalled" aligner:** the bwa read counter goes flat for 30–50 min
routinely — progress lines are emitted sporadically. **Flat counter + full CPU
(~900%) = slow patch, wait.** The actionable failure is flat counter + *low* CPU.
Check with `/proc/<pid>/stat` fields 14+15 over 8 s. Killing a computing job costs
hours; this happened twice and both recovered on their own.

---

## 10. WHAT TO DO NEXT

**Immediate (blocked only on alignment):**
1. Run `s7c_editor.py` on A3A-Y130F-clone2 vs D10A-clone1 (§7). Validate first.
2. Repeat with clone5, clone7 as replicates — an effect in one clone but not the
   others is clone-luck, which is exactly how earlier phases failed.
3. Verify cross-study germline concordance (§11) before trusting a Parent-masked
   result.

**Then, in rough priority:**
4. **Lj-BE as a cross-family control.** Lamprey CDA is a different deaminase
   family. If the hairpin preference is A3A-specific, Lj-BE should NOT show it.
   Strongest available specificity test; samples already queued.
5. **eA3A-RL1** — engineered A3A vs wild-type-lineage A3A-Y130F.
6. **haA3A (Y130G, VA)** — pre-registered expected-null (engineered for
   near-background off-target activity). A null supports the design; a positive
   would need the §5 traps excluded first.
7. **YE1 (rAPOBEC1)** — enzyme-class contrast, same study as the nCas9 calibrators
   so it is internally controlled.

**Do NOT spend effort on:**
- A DNA language model. The ±5 bp saturation result (§4) says there is nothing
  beyond local context to find. Measured, not assumed.
- Cell-state ssDNA proxies (ATAC, R-loop, Repli-seq). Null in every earlier phase,
  including a matched-cell-type spKAS direct ssDNA map (+0.000).

---

## 11. STILL UNTESTED / OPEN

- **Cross-study germline concordance.** The first editor test may use `Parent`
  (PRJNA1042830) as germline mask for P66 (PRJNA1006866) samples. Our 98.76%
  homozygous concordance was measured **within** PRJNA1042830, and 293T sublines
  drift between labs. **Run the same check — Parent vs a P66 control — as soon as
  any P66 sample has counts.** If concordance is low, use P66-background only,
  which means waiting for it.
- **Site-level coverage matching in enrichment.** Chromosome-level showed no gross
  confound (ρ=+0.258, p=0.23) but that is coarse; the confound that mattered in
  earlier phases operated per site.
- **`pos_in_loop` convention** does not match Buisson's published pattern exactly.
  Low priority; stem length is the load-bearing feature.

---

## 12. BUGS FOUND AND FIXED (all one family)

A convention correct locally, then consumed downstream as if universal. **None
crashed. Each would have returned a confident wrong number.**

1. Negatives sampled plus-strand-only → exactly-ceiling **11.000×** enrichment.
   Caught by noticing 11.000 = 1/base_rate.
2. −2 base via a `ref=='C'` branch → YTCA/RTCA and the A3A/A3B donor split
   corrupted. Caught by deriving the genomic background (0.6057) and finding
   donors sitting *below* it.
3. universe `pos` 0-based vs trainset `pos` 1-based, same field name. Caught
   latent, before it ran.
4. `s6_pileup` counting only uppercase `T` → would drop half the reads
   **strand-asymmetrically** (the old Selict 21.8:1 mechanism).
5. `pos1` missing from 22 of 23 universes — the bug-3 fix patched the artifact,
   not the generator.
6. The reference-base guard was **vacuous** — `mpileup -f REF` reports the base
   from the FASTA, not the reads, so it compared a FASTA to itself. An hg38 BAM
   against the hg19 universe gave 0/158,973 mismatches and would have PROCEEDED.
   Replaced with a BAM `@SQ`-vs-`.fai` build check.

**Also corrected: two of my own framing errors.** (a) I applied the ≥3× bar to
enrichment when it belongs to burden — that would likely have called a
textbook-behaving A3A editor a null. (b) I recommended restarting P66-background
at 9 threads to "save 13 h"; threads are zero-sum on a 32-core box, so the real
lever was queue order, not thread count.
