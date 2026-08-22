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
# FIRST RUN -- only D10A-clone1 exists at that point, and P66-background lags,
# so Parent is the germline mask. Extra controls MUST be omitted: s7c intersects
# chromosomes across all named samples, so naming an unaligned control yields an
# empty intersection and the run aborts.
~/miniconda3/envs/apobec/bin/python /mnt/data/a3a/s7c_editor.py \
    P66-A3A-Y130F-clone2  P66-D10A-clone1  ''  Parent

# LATER, once D10A-clone6/clone10 and P66-background have 23/23 counts:
~/miniconda3/envs/apobec/bin/python /mnt/data/a3a/s7c_editor.py \
    P66-A3A-Y130F-clone2  P66-D10A-clone1 \
    P66-D10A-clone6,P66-D10A-clone10  P66-background

# args: <editor> <calibrator> <extra_controls_csv|''> <parent>
# GATE: run s8_xstudy.py Parent <any P66 control> FIRST. If hom concordance
# < 0.90 the Parent mask is invalid for P66 samples -- wait for P66-background.
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

---

# 13. OVERNIGHT RUN 2026-08-20/21 — RESULTS

Written during a gcloud auth outage from the session's own logs. Every number below
was produced and verified on-node before the outage; nothing here is reconstructed.

## 13.1 THE EDITOR TEST RAN. Endpoint B: NULL. Endpoint A: INCONCLUSIVE.

Decision table was fixed in writing at 09:45 UTC, **before any P66 number existed**
(STATUS.md, "PRE-REGISTRATION"). Results read against it, not the other way round.

**STEP A — cross-study germline gate. PASSED.** This was §11's top open question.
```
s8_xstudy Parent vs P66-D10A-clone1, 23 chromosomes
  germline HOM concordance  0.9839   (within-study reference 0.9876)
  germline HET shared       0.4720   (within-study 0.4942)
  LOW-VAF 0.03-0.15 shared  0.0637   (within-study 0.0677 — LOWER = safe direction)
```
≥0.97 → same lineage. **The two labs did not drift their HEK293T sublines.** Parent is
a valid germline mask for P66. The low-VAF check also passes in the safe direction; a
HIGH cross-study low-VAF share would have implied a shared technical artefact.

Context for reading 0.9839: two *same-study, same-lineage* nCas9 clones score only
0.9706 against each other when both are shallow (measured 03:45). Cross-study 0.9839 is
comfortably inside same-lineage territory.

**STEP B — validation. PASSED bit-for-bit.** nCas9-clone2 vs nCas9-clone1:
burden 0.871 / 0.879 / 0.877 / 0.886, cov ratios 1.005 / 1.002 / 0.999 / 0.997.

**STEP C — editor. P66-A3A-Y130F-clone2 and -clone5 vs P66-D10A-clone1, Parent mask.**

ENDPOINT A — BURDEN: **INCONCLUSIVE, and the calibrator is disqualified.**
Only one of four coverage bands was depth-matched (cov ratio 0.989); the other three
were auto-flagged UNMATCHED by the depth audit. In the matched band the ratio was
**0.144** (clone2) and **0.137** (clone5) — the editor carrying ~7× FEWER specific
sites per Mb than its deaminase-free calibrator. Direction inverted, not merely short
of the ≥3× bar.

Cause found:
```
sample                med_cov     alt>=1     alt>=2  alt>=2/Mb  germline VAF>0.9
A3A-Y130F-clone2         27      864,617    184,003     794.5       61,174
A3A-Y130F-clone5         31      141,625*        —          —            —
D10A-clone1              41    6,457,037    339,969    1469.0       60,240
Parent                   35    1,135,489    307,860    1332.6       60,555
nCas9-clone1             26      854,914    256,302    1117.4       60,537
                                                   (*at cov>=20, chr1+2+22)
```
D10A-clone1 has **7.5× more alt≥1 sites than any other sample**, and **94.8% of them sit
at VAF<0.05** (single-read band) against 78–79% for the editor clones. Its VAF≥0.35
share is 1.4% vs 10.6–11.4%. Contamination or a real subclone would populate 0.05–0.35;
that bin holds 0.4%. Germline VAF>0.9 counts are near-identical across all samples, so
this is not lineage or calling drift. **D10A-clone1 has a sequencing-quality defect and
is not usable as a burden calibrator.** At its median depth of 41 the mechanism reaches
alt≥2 as well, which is the direct cause of the inverted ratio.

ENDPOINT B — HAIRPIN (within-sample ratio, depth-immune → the interpretable half):
```
stem   clone2   clone5   pooled   pooled null95   n_hp   calibrator
  6    0.899    0.859    0.879     0.93-1.07      589      1.213
  7    1.034    1.144    1.089     0.88-1.14      250      1.378
  8    1.069    1.323    1.196     0.78-1.22       96      1.532
```
**Neither editor clone reaches the calibrator at any stem length.** Pooled stem-6 is
BELOW its null lower bound — a significant *depletion*, the opposite of a
hairpin-targeting signature.

Clones are NOT distinguishable from each other (Fisher p = 0.59 / 0.45 / 0.31), so the
pooled estimate is the right one. Both editor clones are technically clean and
indistinguishable (79.0% vs 78.5% VAF<0.05; 11.4% vs 11.4% VAF≥0.35).

POWER: the calibrator's value lies OUTSIDE the editor's own null95 at every stem
(1.213>1.11, 1.378>1.19, 1.532>1.34). The editor **could** have detected a
calibrator-sized effect. This is evidence of no calibrator-sized effect, not silence.

## 13.2 **1642× CROSS-CLONE RECURRENCE — and it was carrying the stem-8 signal**

New control, not in any prior phase. Two *independent* clones should share ~2
editor-specific sites. They share **3,000**:
```
jointly eligible   217,045,974
clone2 specific         19,975
clone5 specific         19,847
shared by BOTH           3,000
expected by chance           1.8
OBS / EXP              1642.45x
```
Not lineage germline — VAF median 0.083 (germline sits near 0.5 or 1.0). At coverage ~30
that is alt=2–3 reads at recurrent positions: systematic error-prone loci, or a
sub-clonal population predating clone isolation. **Not editor-attributable either way.**

Splitting clone2's sites by whether clone5 also called them:
```
set        stem   n_hp   n_sites    enr
shared        6     40     3,000   0.795
shared        7     22     3,000   1.275
shared        8      8     3,000   1.326
private       6    262    16,975   0.920
private       7     97    16,975   0.994
private       8     35    16,975   1.026
```
**On clone-private sites the editor is null to three decimals at stems 7 and 8.** The
marginal pooled stem-8 value (1.196, p=0.046) came from the shared fraction. Note the
pooled n_hp at stem 8 is 96 — the A3A-vs-A3B claim in this project died at n=97.

ACTION TAKEN: `s7c_editor.py` gains `A3A_EXCLUDE_RECURRENT=<siblings>`, default OFF,
and every run prints whether the filter was active and how many sites it cut. Patch
written; apply with `fix_s7c_recurrence_v2.py`.

## 13.3 MODEL TRACK — PURITY BEATS VOLUME (matched N, matched block)

```
arm   tcw     N        purity   AUROC             top0.1%  RANDOM   n
v4s  >=0.20   83,999   0.8158   0.5506 ±0.0091    3.933x   0.965x   923
v5   >=0.40   83,999   0.8642   0.5382 ±0.0080    5.017x   1.061x   923
```
Identical N, identical 6-feature block, identical fold structure. **+27.6% tail for
+0.049 purity.** v5's fused block reaches **5.280×**, beating v2's 4.960× with 4.1×
LESS data.

And the converse: **more data HURTS when it is dirty.** v3 (1,000,000 positives,
uncapped, hypermutator-heavy) gives the HIGHEST AUROC (0.6115) and the WORST tail
(**2.904×**). Selecting a dataset by AUROC picks the worst one available.

**Full v5 tail ladder:**
```
gc_only                 1.120x   RANDOM 1.061x   ← inert
sequence_only           3.814x
hairpin_nogc            4.839x
hairpin_only            5.017x
hairpin_nogc+sequence   5.148x
hairpin+sequence        5.280x
```

**ARCHITECTURE — this REFINES §10's "do not build a DNA LM".** Separate the claims:
- **Size**: still argued against, three independent ways — context saturates by ±10 bp
  (tail *decays* 2.351× at ±5 bp → 2.274× at ±40 bp), ~100–400k usable positives not
  2.38M, and v3 showing more dirty data lowers the tail.
- **Structure branch**: now POSITIVELY SUPPORTED. On v5, sequence_only 3.814× →
  hairpin+sequence **5.280×** = **+38% tail for +0.0023 AUROC**. A model selected on
  AUROC would discard the structure branch entirely.

**Corrected recommendation: build the dual encoder, but COMPACT** — short-context
sequence encoder (no Evo/NT-scale), plus the structure branch, trained on a CURATED
high-purity set, **selected on tail enrichment, never AUROC.**

Realized purity WITH the burden cap is **0.8155 / 0.8642**, not the 0.869/0.939 quoted
in the cron prompt — those are UNCAPPED. The cap and the purity filter pull in OPPOSITE
directions (at tcw≥0.40 the cap LOWERS purity 0.8909→0.8642) because high-tcw donors
are hypermutators. Nobody had computed this.

## 13.4 CONFOUNDS TESTED AND CLEARED THIS RUN

- **GC deciles**: hairpin survives 10/10 with in-stratum random baselines.
- **Complexity**: survives 10/10 entropy deciles, and the DIRECTION refutes the
  artefact hypothesis — enrichment is WEAKEST in the lowest-complexity decile (1.071×)
  and strongest in high-entropy deciles (1.598×). Dropping the most repetitive 50% of
  the data moves it 1.2705 → 1.2571.
- **Hotspots**: 99.93% of positive coordinates are singletons; recurrent coords account
  for 0.14% of positives.
- **Hairpin feature verified against hg19**: 65,000 sites, stem-pairing / gc_pairs /
  hp_score / focal-C-in-loop all **100.000%**, including strand=1 sites (the bug-2
  orientation test). Search is strand-symmetric to ≤2%.
- **Coordinates**: ref base at 0-based `pos` is C/G by strand in **100.0000%**; using
  `pos1` as the index gives 0.0000% — the two conventions are maximally distinguishable.
- **−2 YTCA channel**: hairpin survives in BOTH YTCA and RTCA strata, and P(stem≥6) is
  equal across them (2.216% vs 2.224% among positives) — the channels are orthogonal.
- **YTCA background 0.6057 derived twice independently**: 0.6063 (83 donors) and
  0.6056 (133 donors).

## 13.5 BUGS 6 AND 7, AND WHAT THEY COST

**BUG 6 — output filename did not follow the input version.** Every sed-derived variant
(`s3_v3/v4/v4s/v5`, `s5_v4`) rewrote only the INPUT trainset name; the OUTPUT was
hardcoded `s3_results_v2.json`. All arms overwrote each other, and every surviving file
was labelled with the *script's* version rather than the *data's*. Caught with the true
v2 baseline minutes from destruction. FIX: output name is now derived from `TRAINSET`,
so it cannot disagree with the input.

**BUG 7 — `hairpin_only` was not hairpin-only.** `local_gc` sat inside the block. The
positive rate slides 1.210× → 0.789× across its deciles (negatives are not GC-matched),
so every `hairpin_only` number was hairpin PLUS an uncontrolled regional channel. FIX:
decomposed into `hairpin_nogc` / `gc_only` / `hairpin_only`. **Verdict: it changed
nothing** — 4.249× → 4.264× without GC, and `gc_only` alone is 1.135× (1.120× on the
pure set, against a 1.061× random baseline, i.e. inert). The defect was real; it was
not inflating the result.

**Also: a landmine I created and removed.** Fixing a hardcoded `envs/bio/samtools` on
ai-chem2 (which broke s6 there for two driver sweeps) left the two nodes with DIVERGENT
copies of `s6_pileup.py` — and ai-chem has samtools ONLY in `bio`. Syncing ai-chem2's
copy onto ai-chem would have broken the critical path. Unified via `_find_samtools()`;
normalised md5 now identical on both nodes.

**samtools versions differ across nodes** (ai-chem 1.24 / ai-chem2 1.23.1). Bounded:
the mpileup invocation is byte-identical and pins `-q 20 -Q 20 -d 500 --no-BAQ`. BAQ —
the largest version-drift source — is explicitly disabled. Residual: `--ff` is not
pinned; only matters for a cross-node comparison, and none has been run.

## 13.6 OPEN, IN PRIORITY ORDER

1. **P66-D10A-clone6** — clean replacement calibrator. Was 94.3% at 11:43 UTC; counts
   likely completed ~12:40. **This is the single highest-value item**: it makes Endpoint
   A measurable for the first time and gives Endpoint B a second clean reference.
   Its BAM is protected by unit `a3a-keeper2` (hardlinks D10A BAMs on sight, because the
   running driver's `KEEP_PATTERN` predates the calibrator's importance).
2. **P66-D10A-clone10** — second clean calibrator. If clone6 and clone10 agree, the
   baseline is established; if they disagree, D10A clones vary technically and the whole
   calibration approach needs rethinking.
3. **P66-A3A-Y130F-clone7** — third editor clone; moves stem-8 past n=96.
4. **P66-background** — matched germline mask, removes the cross-study Parent step.
5. Apply `fix_s7c_recurrence_v2.py` and rerun both editor clones vs clone6, **with and
   without** the filter, so its effect is visible rather than assumed.
6. `qa_specset_v2.py` replaces `qa_specset.py` (v1 has dead code and a fragile init).

## 13.7 OPERATIONAL LESSONS

- **The percentage counter in the cron prompt was ~2× wrong all night.** bwa logs are
  append-only and were created 2026-08-10; each holds a dead 10-day-old run PLUS the
  live one, concatenated with no banner. Summing the file adds them. Tell: two batch
  sizes in one log (600000 = dead generation, 466668 = live). Replaced by
  `progress.sh`, which reads bwa's own open-FASTQ fd offset — no history, and the two
  mate files agree to ~0.1% as a built-in consistency check.
- **A niced background job is not automatically harmless.** Nice governs CPU shares, not
  cache or memory bandwidth. A 23.5 GB-resident model job with 48 threads roaming 32
  cores halved bwa throughput at only 12% CPU. `taskset` to 4 cores recovered 1.8×.
- **Eight slow-patch episodes**, 11–25 min each, all self-resolving: bwa at full CPU
  consuming no input, because FASTQs are in flowcell order and hard reads arrive in
  contiguous blocks. Never a stall. Confirm with per-process CPU before acting.
- **CHECK THE INSTRUMENT BEFORE THE MACHINE.** Eight of my own measurement errors this
  run (self-matching `pgrep`, misindexed probe, `awk` int32 overflow, grep pattern
  hiding fold lines, a rename loop ignoring its own check, testing impossible input,
  reading the wrong file for a systemd unit twice). The system was right every time.
- **Resolve `systemctl show <unit> -p ExecStart` before reading any script.** I audited
  the wrong file twice; both times the stale sibling was the weaker version. Dead scripts
  are now renamed `SUPERSEDED_*`.
- `feat/ARTIFACTS.md` on both nodes is the authority on what each result file actually
  contains. If a filename and the manifest disagree, trust the manifest.

---

## 14. DeaminaFormer architecture search — 2026-08-21/22 overnight (COMPLETE, decisive)

The question was: does a pretrained DNA foundation model, fused with structure and sequence
features, beat the gradient-boosted baseline at the top 0.1% tail? Ran end to end on a
V100 (`ai-gpu`, us-central1-c, 200 GB scratch PD at `/mnt/a3a`). **Answer: no, and the
reason is now pinned down rather than guessed.**

### 14.1 The decisive table
Same data (`a3a_trainset_v5`, 923,989 rows, base rate 1/11), same 5 folds by **held-out
chromosome**, same pooled-OOF top-K% metric, random baseline beside every number.

|                              | MLP head | GB head |
|------------------------------|---------:|--------:|
| struct only (6 feat)         |  4.881×  | 4.940×  |
| struct + thermo (16 feat)    |  4.655×  | **5.036×** |
| NT-v2 embeddings alone       |  1.250×* | 1.048×  |
| NT-v2 + struct               |  1.488×  | 4.857×  |
| everything                   |  2.333×  | 4.774×  |
| GB one-hot + hairpin (baseline) | — | 5.280× |

random 0.917 · arithmetic ceiling 11.000 · n=924 at top 0.1%

\* Both ladders ran an ntv2-only block under an MLP head: the flat-MLP ladder scored
1.155× and the CNN ladder's control scored 1.250×. The table quotes the latter. Both sit
against a random baseline of 0.917, so the conclusion is identical either way.

### 14.2 Three conclusions, and the control that separates them
1. **The embeddings are EMPTY for this task, not merely diluted.** GB subsamples features
   and can ignore noise, yet 128 PCs of NT-v2 score 1.048× against a random baseline of
   0.917. Two heads, two architectures, same answer.
2. **The MLP's collapse was the HEAD.** Adding embeddings to structure cost the MLP almost
   everything (4.881 → 1.488) and cost GB essentially nothing (4.940 → 4.857). The
   encoder-vs-head control was queued **before** any of these numbers were seen — without
   it, "embeddings hurt" and "embeddings are empty" are indistinguishable.
3. **Nothing beats hairpin geometry.** A Conv1d motif scanner over ±200 bp (2.5× the
   baseline's window, the right inductive bias) adds nothing to structure (4.810 vs 4.881)
   and reaches 2.595× alone.

**AUROC was inverted against the tail at every step.** Among CNN blocks `cnn_only` has the
*best* AUROC (0.5671) and half the tail of `struct_only` (worst AUROC, 0.5388). Selecting on
AUROC would have chosen the weakest model at every decision.

### 14.3 The recommended model
**16 features (6 hairpin geometry + 10 thermodynamic), gradient-boosted.**
- top 0.1% enrichment **5.036×** (random 0.917)
- **ECE 0.00125** (spec < 0.05) — and calibrated at the operating point, not just on
  average: at the top 0.1% it predicts **0.4822**, observes **0.4782**
- category-agnostic: 4.4–5.4× in CDS, intron/UTR and intergenic alike
- auditable, cheap, mechanistically interpretable — a better outcome for a regulatory gate
  than a foundation model would have been

### 14.4 Supporting results
- **A model is required; a threshold is not enough.** Best single hand rule (no training,
  no folds, which *favours* it) is hp_score at 2.988×. stem≥8 covers 2,229 sites at a
  1.860× positive-rate ratio; the top 0.1% needs 924, so every threshold must cut inside a
  tied block. Ranking within ties is what the model supplies.
- **Site identity is not learnable.** Disjoint donor halves (10 vs 11 donors) share 10 sites
  where chance predicts 4.6–17.6 depending on the assumed TCW universe — **at chance**.
  99.974% of sites appear in exactly one donor. The learnable target is a per-site
  *propensity*, never *which* site is hit.
- **The thermodynamic block** (DNA partition function, Mathews2004 **DNA** parameters —
  RNA params give MFE −21.1 vs DNA −8.9, a 2.4× error) scores 2.476× alone and adds ~2%
  to struct under GB. `p_unpaired` largely re-encodes stem length.
- **Confounds cleared on v5**: GC 10/10 deciles; complexity 10/10 on two independent
  measures (dinucleotide entropy and distinct-4-mer count), both *weakest* at lowest
  complexity, refuting the repeat-artefact reading; coverage/mappability — MH-adjusted
  1.323 vs crude 1.349 at stem≥8, and the model tail survives coverage-matched negatives
  (5.280 → 5.096, inside the noise floor).
- **Hairpin features are class-consistent.** An independently written inverted-repeat
  scanner relates to the stored `stem` identically in both classes (Δcorrelation 0.0096,
  Δmean-offset 0.0314) — this is the check the whole result rests on, because a
  class-dependent feature would have fooled both learners identically.

### 14.5 Empirical noise floor — read every comparison against it
The built-in random baseline should be exactly 1.000×. Across 18 blocks it is
**0.973 ± 0.062** (range 0.870–1.085). **Differences below ~12% at top 0.1% are not
resolvable.** This retired one earlier claim ("v5 fused *beats* v2" → *matches* v2) and it
covers the 5.036 vs 5.280 gap.

### 14.6 Bugs found and fixed in this block of work
- **BUG 6, three instances**: sed-copied scripts writing their parent's output name. Found
  in `.json` outputs (4 scripts) and then again in `.npy` OOF outputs (8 scripts) after I
  had declared the family closed twice. *"I fixed that bug family" is a claim about a
  SEARCH*, and mine was incomplete twice.
- **N-handling, three instances in one night**: `np.clip(win,0,3)` silently recoding N→T;
  NT-v2's 6-mer tokenizer falling back to single-character tokens on N (a 1,025-N window
  becomes 1,025 tokens, padding the whole batch → the CUDA OOM at 96%); and the same
  ambiguity in the one-hot path. All now encode N as *absent*, never guessed.
- **Head-slice on block-ordered files, three instances in my own QA code.** These files are
  written positives-first, so `[:200000]` samples 45% positives instead of 9%. Fixed by
  random sampling **plus a printed base-rate assertion** — writing the lesson down did not
  stop me repeating it; a check did.
- **Host OOM killed the MLP ladder** on its 1,620-dim final block after I had fixed that
  exact fragility in the sibling script two ticks earlier. Results recovered by parsing the
  log; generators now dump after **every** block.
- **Two gate races**: `chain3/4` polled for a process that had not started yet; the
  replacement waited on an artifact that **already existed from an earlier run**. Both would
  have run two ladders on one GPU. Fixed with a timestamp gate (`touch` a stamp, require the
  artifact to be *newer*).
- **`pkill -f <pattern>` matched my own SSH command line twice**, killing the session.
  Bracket form (`"[c]hain.sh"`) or kill-by-PID.

### 14.7 What is NOT done
- **Tier A/B stratified recall — the spec's headline metric — is blocked on data**, not
  analysis: COSMIC Tier 1 (581), ClinGen HI Level 3 (~340), DepMap common essentials
  (~1,800). None on either node. Only the C/D split (refGene) could be built.
- **Calibration is to the sampled 1-in-11 base rate**, which is a sampling artefact, not
  genome-wide prevalence. Any deployed threshold must be re-based before a risk number is
  quoted.
- **v5 rests on 21 donors**, not the 133 available at v2 — purity cost donor diversity, and
  that bounds generalisation.
- Fine-tuning the foundation model was **pre-registered as not warranted** if frozen
  embeddings failed to beat GB. They scored 1.048×. Not attempted, deliberately.

### 14.8 Where the artifacts are
- `results/deaminaformer_ablation.json` — flat-MLP ladder, 11 blocks (`complete: false`,
  names the OOM-killed block)
- `results/deaminaformer_v2.json` — CNN ladder, 9 blocks
- `results/gb_on_embeddings.json` — the encoder-vs-head control
- `scripts/` — every generator, including the QA scripts that produced the confound clearances
- On `ai-gpu`'s scratch PD (persists across the stop): `emb_ntv2_v5.npz` (1.9 GB),
  `emb_hyena_v5.npz` (0.5 GB), `struct_thermo_v5.npz`, `a3a_trainset_v5_win1k.npz`
- `ai-gpu` is **STOPPED**; no foreign process was running on it.

---

## 15. The editor track has a hard measurement bound — 2026-08-22

Five arms across two studies all landed at or below the deaminase-free calibrator. The reason
turns out not to be biology, and it is worth stating before anyone runs more clones the same way.

### 15.1 Endpoint B is a subclonal measurement by necessity, not by choice
The Parent-masked, `alt>=2` specific-site sets are **91–94% below VAF 0.15**. And hairpin
enrichment **reverses across VAF strata**: at VAF<0.05 it is 1.04–1.15×, at VAF≥0.15 it is
**0.545–0.660×** (depleted). Two studies, two deaminase families, same pattern.

### 15.2 The clonal test is not weak — it is impossible
Background P(stem≥8) on the Parent-silent set is 0.00201, so among a sample's clonal variants
the hairpin count **expected by chance** is n × 0.002:

| sample | all spec | ≥0.15 | ≥0.35 | expected n_hp |
|---|---:|---:|---:|---:|
| A3A-Y130F-clone2 | 20,010 | 2,710 | 601 | 1.2 |
| A3A-Y130F-clone5 | 19,899 | 2,477 | 462 | 0.9 |
| nCas9-clone1 | 90,844 | 7,000 | 554 | 1.1 |
| nCas9-clone2 | 80,053 | 6,843 | 639 | 1.3 |
| D10A-clone1 | 172,809 | 5,010 | 337 | 0.7 |
| Lj-BE-clone3 | 1,301,573 | 28,405 | 542 | 1.1 |
| Lj-BE-clone5 | 203,390 | 4,020 | 305 | 0.6 |
| eA3A-RL1-clone1/2/5 | 186k–498k | 3.3k–8.3k | 352–396 | 0.7–0.8 |

**One expected site per sample**, controls included. A test whose null expectation is a single
site cannot resolve a 1.3× effect, let alone a 3× one. This is a bound, not a null.

**For a DNA safety gate this is the crux**: the clinically relevant off-targets are the
clonally fixed ones — permanent and heritable. The assay as constituted cannot address them.

### 15.3 What it would take
| target expected n_hp | clonal variants needed per arm | |
|---:|---:|---|
| 10 | 4,975 | barely |
| 30 | 14,925 | thin but usable |
| 100 | 49,751 | solid |

Current: **337–639 per sample**. So ~25× more for a thin test, ~150× for a solid one.
**Depth does not buy this** — a variant is clonal or it is not. It needs **25–50 clones per
arm instead of 2–3**, or a different assay design.

### 15.4 Method changes this forced, all validated
- **`alt>=2` is a depth-dependent VAF cut** (VAF≥0.05 at 40×, ≥0.087 at 23×), so arms of
  unequal depth were never comparable. A **VAF floor failed its own validation** — two
  deaminase-free clones diverged 1.438 vs 1.159. **Matched coverage bands passed** it
  (1.118/1.081, 1.383/1.380, 1.583/1.640) and were confirmed against an independent
  re-implementation, 6/6 values to three decimals. Banded endpoint B is the method of record.
- **The noise-floor criterion must be applied WITHIN study.** Normalised per × of coverage,
  node A sits at 33k–39k and node B at 143k–466k — a 4–14× gap after normalisation, so the
  studies genuinely differ in error rate. D10A-clone1 (205,900/×) sits *inside* node B's band:
  it looks like a sample from the noisier study placed in the quieter cohort. That is a
  sharper account of its defect than "7.5× outlier", and applying node A's threshold to node B
  wrongly disqualified two sound Lj-BE clones for about twenty minutes.
- **Calibrator qualification runs before any editor number** (`qualify_calibrator.py`,
  wired as Step A0). Validated to reject D10A-clone1 and pass all five sound samples — the
  first version rejected a *good* calibrator and was caught by testing it against known-good
  samples, not just the one it was meant to reject.

### 15.5 Corrections made to the record
- **"The shared fraction carried the entire stem-8 signal"** (A3A-Y130F) rests on **n_hp = 8**;
  Fisher p = 0.520. The downgrade still stands, but on the other leg: 3,000 sites shared
  between two independent clones against 1.8 expected (1642×) means those sites are not
  clone-private mutation, whatever their hairpin content.
- **The stem-6 between-arm difference** (A3A-Y130F 0.877 vs eA3A 1.061, p = 2.6×10⁻⁶) reverses
  under VAF control. It measured VAF composition. Retracted.
- **"94.8% of D10A-clone1 sites below VAF 0.05"** is 87.2% under an explicit `alt≥1 & cov≥8`
  definition, and 31.2% on the set the editor test actually uses. The figure moves with its
  eligibility filter and was being quoted without one.
- The **6.4× site-count gap** between Lj-BE clone3 and clone5 is **subclonal load, not noise**:
  the excess peaks at VAF 0.10–0.15 (10.8×) and vanishes at clonal VAF (1.83×).

### 15.6 A second compositional confound on endpoint B: GC — found 2026-08-22 05:35

Coverage banding fixes **depth**. It leaves **GC composition** entirely untouched, and that is a
confound of the same size and shape. Check 4 had been run all night on the PCAWG model and
**never once on the editor arms**.

The hairpin background is strongly GC-dependent — inverted repeats are far commoner in AT-rich
sequence:

| GC decile | p_bg(stem≥6) |
|---|---:|
| 0.012–0.284 | 0.03564 |
| 0.543–0.951 | 0.01248 |

A **2.9× span**. So any shift in the GC composition of the called sites moves the pooled
enrichment directly. Stratifying A3A-Y130F-clone2 against nCas9-clone1:

| | crude OR | GC-adjusted MH | shift |
|---|---:|---:|---:|
| editor | 0.895 | 0.760 | **−15.0%** |
| calibrator | 1.053 | 1.041 | −1.2% |

**The confound is asymmetric** — 15% on the editor, 1.2% on the calibrator — and 15% exceeds the
~12% empirical noise floor. An asymmetric compositional effect of that size is exactly what
manufactures a spurious between-arm difference, and one such difference has already been
retracted tonight for the analogous reason on VAF.

**Endpoint B requires GC stratification in addition to coverage banding.** Neither substitutes
for the other. `ops/auto_advance_gc.sh` runs it on the clone6 pairs plus the deaminase-free
pair as a null control for the estimator itself.

Direction here is unaffected (adjustment makes the editor *more* depleted), but the method
point holds regardless of which way it happened to move.

### 15.7 An instrument note worth carrying forward
`ops/progress.sh` accumulated **nine** defects in one night, every one found because a number
it printed did not survive a second look — a hardcoded batch size from the wrong node, a
read counter that double-counts re-aligned samples, a merge row timing the alignment instead
of the merge, and four separate failures of a recent-rate ETA.

The ETA was eventually **deleted rather than fixed a fourth time**. Its scoreboard was four
false alarms and zero true positives, and the cause was structural: batches take 8–25 min and
the column was sampled every ~10 min, so a two-point rate is dominated by where batch
boundaries fall. It now reports **minutes since last output** — a direct measurement that
cannot produce a 103.9h reading, and the signal every intervention decision actually used.

Eight suspected stalls were investigated tonight and **all eight were slow patches that
self-resolved** (11–25 min, hard reads arriving in contiguous flowcell-ordered blocks). The
cost of checking was ~2 min each; the cost of acting on one would have been ~16 h, since the
worker writes a `FAILED` flag on a killed pipeline that **blocks automatic re-queueing**.

### 15.8 A queue-level gap: the haA3A arm has no deaminase-free control

All thirteen `SRR257250xx` samples are **one study** (PRJNA1006866): node A's A3A-Y130F ×3,
D10A ×3 and background, plus node B's eA3A-RL1 ×3 and Lj-BE ×3. So **D10A-clone6/clone10 are
the within-study calibrators for the node-B arms too** — the analyses run tonight used
eA3A-clone1 as the comparator for Lj-BE, which is editor-vs-editor and was labelled "calib".
Informative for the cross-family question, but not a calibrator comparison.

The six queued `SRR26881532`-series samples (Y130G ×2, VA ×2, YE1 ×2) are a **third study with
no deaminase-free control in either queue**. Since noise floors differ 4–14× between studies
*after* depth normalisation, borrowing D10A or nCas9 would reproduce the cross-study confound
that made the 05:20 preview uninterpretable.

**So the standing haA3A pre-registration cannot be tested as specified**, and ~18 h of
alignment is queued behind samples whose primary endpoint has no valid comparator. Either
locate a deaminase-free control from that BioProject, or run the arm descriptive-only with the
limitation stated up front.

---

## 16. THE TRANSFER TEST — the project's central claim, tested 2026-08-22 06:40. It fails.

Does a PCAWG-trained classifier rank real base-editor off-target sites above matched
background in HEK293T? Until now the project had (a) the model predicting held-out PCAWG
chromosomes at 5.036×, which is tumour mutations predicting tumour mutations, and (b) the
hairpin *feature* measured separately on both datasets. **Neither is the model transferring,
and the model had never been applied to a single HEK293T site.**

**Design.** Train on PCAWG v5 (the 6 structural features that won the architecture search),
score every jointly-eligible TCW site in HEK293T (Parent-silent, cov≥8 in both samples), and
ask whether editor-specific sites are enriched in the model's top K%. **The calibrator is the
control**: if both arms enrich equally, the model predicts *where variants get called* — the
coverage, mappability and GC artefacts characterised in §15 — not anything about the editor.

### Result, both editor clones, ~216 M sites scored each

| top | editor | n_hit | calibrator | n_hit | ed − cal |
|---|---:|---:|---:|---:|---:|
| **clone2 vs nCas9-clone1** | | | | | |
| 0.1% | 0.806× | 16 | 1.289× | 117 | −0.483 |
| 1.0% | 0.721× | 143 | 1.089× | 988 | −0.368 |
| 5.0% | 0.954× | 946 | 0.962× | 4,365 | −0.008 |
| **clone5 vs nCas9-clone2** | | | | | |
| 0.1% | 0.861× | 17 | 1.063× | 85 | −0.201 |
| 1.0% | 0.628× | 124 | 1.003× | 802 | −0.374 |
| 5.0% | 1.001× | 988 | 0.997× | 3,989 | +0.004 |

**Positive control, same fitted model on PCAWG in-sample: 5.369× at top 0.1%.** The model
works. What fails is the transfer.

### What this establishes
1. **The PCAWG-trained model does not predict base-editor off-targets.** At the
   well-powered depth (top 5%, n_hit ≈ 950–990 per arm) both clones sit at **0.954× and
   1.001×** — flat. At the sharp tail it is *depleted*, 0.63–0.86×.
2. **The editor−calibrator gap is negative at every depth but the last.** The model ranks the
   *deaminase-free control's* sites higher than the editor's — the opposite of the hypothesis.
3. **The mild calibrator enrichment (1.0–1.29×) is the calling-bias signal**, exactly what the
   control was included to expose.
4. **It replicates across two independent clones** with independent calibrators.

### Limits, stated plainly
- The editor-specific set is **91–94% subclonal** (§15.1), so this tests prediction of
  subclonal calls. The clonal test remains impossible at this clone count (§15.2).
- Both calibrators are cross-study (PRJNA1042830 vs PRJNA1006866). The within-study rerun
  against D10A-clone6 is queued and should be done before this is written up.
- n_hit at top 0.1% is 16–17; the top-5% row is the one to read.

### 16.1 The transfer result audited against its best alternative explanation — it survives

A strong negative deserves a positive's scrutiny. The worry: the model ranks hairpin-rich
sites high, hairpins favour AT-rich sequence, AT-rich sequence maps poorly — so the top-ranked
sites might simply be where variants are hardest to call.

- **In-distribution**: the scored universe is TCA 0.472 / TCT 0.528 and nothing else, matching
  the training set exactly.
- **No coverage confound, opposite direction**: the model's top 0.1% has mean coverage **26.66**
  against 25.98 overall. Calling power is not what excludes them.
- **Survives coverage banding**: editor 0.296× / 0.665× / 0.768× / 0.857× across the four bands,
  calibrator 0.788× / 1.138× / 1.053× / 1.217×, gap negative in every band. Best-powered bands
  (n_hit 48, 56) give 0.665× and 0.768×.

**And the audit produced the mechanism.** The model's top 0.1% is **GC-rich (0.4562)** while the
editor's called sites are **AT-shifted (0.3537–0.3784)**. What PCAWG taught the model about
APOBEC mutagenesis in tumours is a sequence preference that does not describe where this editor
deposits damage in HEK293T. The preference is non-monotone in score (top 5% is AT-rich at
0.3554), so "the model likes GC" is too simple a summary — but the mismatch at the operating
point is the concrete reason the transfer fails.

## 17. Node-B Lj-BE arm — both endpoints null (2026-08-22 06:49–06:55)

Three Lj-BE clones (integrity 0 defects, strand 1.003–1.008) completed and both armed
analyses fired.

**Three-way cross-clone sharing.** ALL-3 class 149,590 observed vs 5.63 expected (26,574×),
2-of-3 63,598 vs 7,408 (8.6×). No Parent mask by design, so ALL-3 is germline-dominated —
consistent with its hairpin enrichment sitting at background (1.109× at stem≥8).

Hairpin enrichment by sharing class, stem≥8 (background p=0.00210):

| class | n | n_hp | enrichment |
|---|---:|---:|---:|
| private | 1,854,295 | 5,378 | 1.379× |
| 2 of 3 | 63,598 | 261 | **1.951×** |
| ALL 3 | 149,590 | 349 | 1.109× |

Read against the **~1.4× deaminase-free floor, not 1.0**: the clone-private class — the one
that would carry real editor mutation — is at or below the floor, and the only elevated number
is the recurrence-prone 2-of-3 class. **Same shape as the A3A-Y130F arm the site-recurrence
control downgraded.**

**GC-stratified banded endpoint B** (vs eA3A-RL1-clone1):

| | crude OR | GC-adjusted MH | shift |
|---|---:|---:|---:|
| clone5 editor | 1.049 | 0.976 | −7.0% |
| clone5 calib | 1.070 | 0.940 | −12.1% |
| clone12 editor | 1.193 | 1.057 | −11.4% |
| clone12 calib | 1.069 | 0.939 | −12.2% |

Both columns move together, which is the signature of a shared GC artefact rather than an
editor difference. Adjusted values sit far below the floor.

**Caveat, flagged before the run and held to:** eA3A-RL1 is an *editor*, not a deaminase-free
control. The "calib" column is editor-vs-editor. This is a cross-family contrast and **no
editor claim rests on it**. The within-family null control (eA3A-clone2 vs clone1) gives the
estimator's own floor and was still running at time of writing.

### 17.1 Methodological result worth carrying forward
GC adjustment moves **every** arm by 7–12%. Crude endpoint-B numbers across this project carry
a GC-composition inflation of about that size, so **any crude OR near 1.1 is consistent with
zero** once GC is controlled.

### 17.2 The within-family null control — it sets the bar, and the bar is clone-luck

eA3A-RL1-clone2 vs eA3A-RL1-clone1: two clones of the *same* editor, so any difference is
clone-luck by construction. Three clones scored against the identical reference sample with
the identical estimator:

| arm | GC-adj MH | n_hp | 1 SE |
|---|---:|---:|---:|
| eA3A-RL1-clone2 — **same family as reference** | 0.947 | 3,273 | 2.5% |
| Lj-BE-clone5 — cross-family | 0.976 | 3,579 | 2.4% |
| Lj-BE-clone12 — cross-family | 1.057 | 6,874 | 2.1% |

**No family effect.** If Lj-BE differed from eA3A as a family, both Lj-BE clones would sit
above the same-family control by a similar margin. Clone5 is +0.029 (inside noise); clone12 is
+0.110. The two Lj-BE clones disagree with *each other* by 0.081 — more than clone5 differs
from the null control. **Lj-BE arm closed, null, on a measured floor rather than an assumed
one.**

**The spread is 5.2× the counting error** (0.110 against a 1-SE counting error of 0.021), so it
is not Poisson and **deeper sequencing will not shrink it**. Clone-to-clone variance sets the
resolution of this endpoint, not read depth — the same wall the editor track hit from the other
direction, where n_hp ≈ 1 per clone forces 25–50 clones per arm.

**Estimator reproducibility — an internal check that passes.** The calibrator column is the
same sample in all three runs and returned crude OR 1.070 / 1.069 / 1.070, MH 0.940 / 0.939 /
0.940. Three independent invocations agreeing to three decimals: the estimator is stable, so
the spread among editor arms is real clone variation, not run-to-run noise.

### 17.3 Pre-registered bar for the critical-path run
Stamped into `logs/auto_advance.log` at 07:45 UTC, **before the driver produced any number**:

> A3A-Y130F-clone2 vs D10A-clone6 must exceed the calibrator by **more than 0.11 in GC-adjusted
> MH OR** — the measured clone-luck spread — to count as an editor effect. Inside 0.11 will be
> reported as null.

The bar is empirical (three clones on one reference) and stricter than "above 1.0".

## 18. Coverage-matched negatives — the last Phase-1–6 confound, controlled; signal survives

Negatives were matched on donor+trinuc+strand but **not** on local coverage/mappability.
`a3a_trainset_v5cov.npz` adds coverage decile to the *same* matching key (preserving trinuc and
strand rather than replacing them) and drops the ~29k negatives that cannot be matched:
923,989 → 895,054 rows, same 83,999 positives, base rate 0.0909 → 0.0938.

Held-out **chromosome** folds; random baseline beside every number:

| block | AUROC | top 0.1% | RANDOM | ratio |
|---|---:|---:|---:|---:|
| hairpin_nogc | 0.5266 | 4.905× | 1.131 | 4.34 |
| gc_only | 0.5187 | **1.000×** | 0.857 | 1.17 |
| hairpin_only | 0.5375 | 4.905× | 1.024 | 4.79 |
| sequence_only | 0.6226 | 4.000× | 1.095 | 3.65 |
| hairpin_nogc+sequence | 0.6232 | **5.119×** | 1.012 | 5.06 |
| hairpin+sequence | 0.6244 | 5.096× | 0.929 | 5.49 |

- **Enrichment survives**: 5.119× against the pre-matching 4.960× — slightly higher, not lower.
- **AUROC survives**: 0.6244 vs 0.6205.
- **`gc_only` gives exactly 1.000×** — GC alone has zero tail concentration, a clean internal
  control and the direct answer to "is the tail just GC?".
- **Ceiling recomputed** as the build script demands: base rate 0.0938 → ceiling **10.656×** (exact 1/mean; a 10.661 computed from
  the rounded base rate appeared in an earlier draft and is superseded),
  not 11.000. The best block is 48.0% of ceiling, so not leakage.
- **Random baselines span 0.857–1.131** at n=895 (±13%), matching the recorded ~12% noise
  floor. Differences below that are not differences.
- **AUROC/tail dissociation again**: `hairpin_only` near-chance AUROC 0.5375 but 4.905× tail;
  `sequence_only` better AUROC 0.6226 but worse tail 4.000×.

**Effect on §16.** The transfer test trained on v5 (unmatched). This shows what v5 learns is not
a coverage artefact, so the 5.369× in-sample positive control rests on controlled ground —
"the model works; the transfer does not" is better supported after this than before.

## 19. Purity vs volume — purity wins by 4.4×, and the third arm is correctly refused

Raising `tcw_frac` shrinks N *and* shifts the donor burden distribution, and burden is an
established effect modifier here (1.821× low-burden → 1.151× hypermutator). A naive purity
sweep would confound purity with burden and with N. The ablation holds the burden cap at
≤10k for every arm and N-matches the comparison:

| arm | donors | top 0.1% | RANDOM | AUROC | n |
|---|---:|---:|---:|---:|---:|
| v4 — tcw≥0.20, N=207,106 | — | 4.554× | 0.966 | 0.6082 | 2,278 |
| v4s — tcw≥0.20, N=83,999 | 83 | 4.338× | 0.870 | 0.6043 | 923 |
| v5 — tcw≥0.40, N=83,999 | 21 | **5.280×** | 1.085 | 0.6246 | 923 |
| v5cov — tcw≥0.40 + coverage-matched | 21 | 5.096× | 0.929 | 0.6244 | 895 |

```
PURITY at matched N       (v5 - v4s) = +0.942x
VOLUME at matched purity  (v4 - v4s) = +0.216x
cost of coverage matching            = -0.184x
random-baseline spread across arms   =  0.215   <- the floor these must beat
```

**Purity is worth 4.4× what volume is worth**, and achieves it with 21 donors instead of 83.
Volume's +0.216 sits exactly at the noise floor and is not a demonstrated effect. Coverage
matching costs −0.184, also inside the floor — effectively free.

### 19.1 The tcw≥0.60 arm is refused, not missing
With the burden cap applied, only **4 donors / 20,979 positives** survive at tcw≥0.60. Counts
quoted without the cap (~152k) don't reflect a runnable arm, and the cap is not optional
because burden is the modifier the design controls for. Running it would produce a result about
four people and label it purity.

### 19.2 Standing caveat on v5
v5 — the trainset behind the transfer test, the architecture recommendation and the 5.280×
baseline — rests on **21 donors**. The effect strengthened as the pool shrank 83 → 21, the
opposite direction from the A3A-vs-A3B claim that died going 53 → 97. Reassuring but not
sufficient: at 21 donors a **leave-one-donor-out spread is the honest error bar** on 5.280×,
and `qa_perfold.py` should be pointed at v5 by donor, not only by chromosome.

## 20. §19.2 resolved — 5.148× is donor-robust (21/21 leave-one-donor-out)

The open caveat: v5 carries the transfer test's positive control, the architecture
recommendation and the 5.280× baseline, and rests on 21 donors. Every error bar quoted from it
was a held-out-**chromosome** bar — spatial stability, saying nothing about whether one
person's mutations carry the tail.

**Read the value carefully — 5.148× appears twice in this document meaning two different
things.** §14's architecture table lists 5.148× as the *85-feature* `hairpin_nogc+sequence`
block. The jackknife below reports 5.148× as the *86-feature* `hairpin+sequence` full-data
reference. These are a numerical coincidence, not a mislabelling: the jackknife script builds
all 6 hairpin features plus 80 sequence features and its log prints `X (923989, 86)`, which
settles which block it ran. At n=923 the enrichment is quantised in ~0.0119 steps, so
collisions are possible; the 86-feature block differing from §14's 5.280× by 11 positives is
ordinary run-to-run variation from early stopping's internal validation split.

```
full-data value              5.148x
jackknife mean               5.110x
jackknife sd                 0.068
min / max                    4.960 / 5.220
spread (max-min)             0.260
JACKKNIFE SE                 0.295     ->  5.148x +/- 0.30
~12% noise floor             0.618
spread / floor               0.42x     ->  DONOR-ROBUST
```

Dropping any single donor moves the tail by **less than half the noise floor**. The most
influential is DO46330 (−0.188, 67,078 rows); the *largest* donor by rows, DO218176 at 91,960,
lands at 5.090× — essentially the mean. **Donor size and donor influence are unrelated**,
which is what you want: the tail is not carried by one person's mutation load.

**Sharpest form.** The spread of the 21 estimates (0.260) is *smaller* than the spread of the
21 random nulls measured alongside them (0.354, range 0.874–1.228). Donor identity is not a
material source of uncertainty here.

### 20.1 The caveat that remains
This is a **jackknife**: it answers "does any single donor carry the result" — no. It does not
answer "would 21 *different* donors give 5.148×". That needs a fresh cohort, and §19 shows the
pool cannot simply be widened without giving up the purity that is worth 4.4× the volume.
**Robust to donor removal; untested against donor replacement.**

### 20.2 Generator fix
The summary block used `vals.ptp()`, removed in NumPy 2.x. It crashed *after* all 21 replicate
lines had printed, so no data was lost and the summary was recomputed from the log. Fixed in
`scripts/qa_donor_jackknife.py` (`np.ptp(vals)`); the copy on ai-chem2 still has the crashing
version and needs the fixed file pushed when access returns.

## 21. The critical-path arm — A3A-Y130F is NULL against a qualified calibrator

**Harness validated first** on two qualified deaminase-free clones: burden 0.871/0.879/0.877/
0.886 (cov ratios 0.997–1.005), matching the pre-registered 0.88 clone floor; endpoint B editor
at or below calibrator at stem 7 and 8. A null below is not a broken pipeline.

**Calibrator substitution, recorded before the run.** Both D10A clones fail gate v3 (5.15× and
2.96× above the depth trend, limit 2.0×). nCas9-clone1 replaces them, with the cross-study
objection measured rather than waved away: hom concordance Parent↔A3A-Y130F-clone2 **0.9871**
against a **within-study reference of 0.9876**; depth 23.60 vs 25.02; depth-adjusted noise floor
1.07× vs 0.99×. Better matched on every axis than D10A ever was.

**Endpoint A — burden (≥3× bar applies):**

| band | editor/Mb | calib/Mb | ratio | cov ratio | n_ed |
|---|---:|---:|---:|---:|---:|
| (25,35) | 79.8 | 490.5 | 0.163 | 0.998 | 7,318 |
| (35,60) | 121.0 | 734.3 | 0.165 | 0.991 | 2,920 |

The two shallower bands self-flagged UNMATCHED DEPTH at cov ratio 1.029 and the harness
excluded them — it polices itself. **This is not "the editor does less damage":** nCas9-clone1
is a high-burden clone (173.9–734.3/Mb) and A3A-Y130F-clone2 a low-burden one (79.0–121.0/Mb).
Endpoint A shows **clone burden dominating**, the documented nCas9 clonal-somatic confound.

**Endpoint B — banded Mantel-Haenszel, the method of record:**

```
stem 6    editor 0.920    calibrator 1.080    ed - cal = -0.160
stem 7    editor 1.080    calibrator 1.381    ed - cal = -0.301
```

The pre-registered bar (stamped into the driver log at 07:45, hours before clone6 finished)
required the editor to **exceed** the calibrator by **more than +0.11** in GC-adjusted MH OR.
Observed −0.160 and −0.301 — **on the wrong side of the bar.** Unbanded agrees: editor
0.899/1.034/1.069 vs calibrator 1.051/1.377/1.625, p_ed 0.974/0.385/0.358.

**Verdict: NULL.** The test was exploratory and returns null; the framing is not being revised
after the fact. n_hp_ed 302/119/43 at stems 6/7/8 unbanded; banded cells 4–60, and the stem-8
banded cells are too thin to pool and are not quoted.

### 21.1 Caveat added on audit — endpoint B may be burden-affected too

§21 said endpoint A is burden-confounded and read endpoint B as clean. I tested that split.

nCas9-clone1 carries **6.1× the burden** of A3A-Y130F-clone2. Two nCas9 clones differing 12.3%
in burden give near-identical enrichment, which first looked like proof that burden doesn't
drive enrichment — **that over-claimed from n=2**. The implied elasticities are −0.361 / +0.270
/ +0.243 at stems 6 / 7 / 8: **the sign isn't even consistent**, so with two clones the
relationship isn't estimable. Extrapolating each stem's elasticity to the editor's 0.163 burden
ratio predicts 2.023 / 0.843 / **1.046** against observed 0.899 / 1.034 / **1.069** — at stem 8
a pure burden effect predicts almost exactly what was seen.

**This does not overturn the null** — the pre-registered bar was missed in the wrong direction
either way. It does correct the framing: the banded MH endpoint is **depth-matched but not
burden-matched**, and endpoint B may be burden-affected too. Settling it needs deaminase-free
clones spanning a burden range or a burden-matched comparison; neither exists, and D10A-clone10
is same-study/same-protocol as the two already disqualified.

### 21.2 Two-clone replication, and a measured true-null floor on node A

Both A3A-Y130F arms in `auto_advance_gc.log` had **crashed** at 09:28 on
`counts_P66-D10A-clone6_chr1.npz` — the same `.partial.npz` mis-naming that stopped the
critical-path driver. That bug cost two analyses, not one. Re-run after repair, against
nCas9-clone1 (D10A-clone6 is disqualified):

| arm | ed MH | cal MH | ed−cal | ed shift | cal shift | ed sites |
|---|---:|---:|---:|---:|---:|---:|
| **NULL** nCas9-clone2 vs nCas9-clone1 | 1.077 | 1.042 | **+0.035** | −2.4% | −1.2% | 79,654 |
| A3A-Y130F-clone2 vs nCas9-clone1 | 0.760 | 1.041 | **−0.281** | −15.0% | −1.2% | 19,841 |
| A3A-Y130F-clone5 vs nCas9-clone1 | 0.739 | 1.041 | **−0.302** | −14.5% | −1.2% | 19,694 |

**Two independent editor clones agree to 0.021** — the first two-clone replication this arm
has had, and not clone luck. Both sit far on the wrong side of the pre-registered bar (>+0.11)
and far below a **measured** true-null floor of +0.035.

**But the script's own diagnostic fires.** GC adjustment moves the editor column −15.0%/−14.5%
and the calibrator only −1.2%: the arms are not compositionally comparable. Note the direction
— adjustment makes the editor *more* depleted, so GC was **masking** part of the depletion, not
manufacturing it. The editor arms also carry **4.6× fewer** specific sites (19.8k vs 90.8k),
the burden asymmetry from §21.1.

**Where it leaves the arm:** the null is replicated and sits outside a measured floor. What is
*not* established is attribution — burden (4.6×) and GC composition (−15% vs −1.2%) both differ,
and this estimator cannot separate either from a deaminase effect.

## 22. Every PRJNA1006866 control fails the noise floor; the editors don't

P66-background reached 23/23 — the same-study control the A3A-Y130F arm was missing. It is
**disqualified at 9.71×** above the peer depth trend (16,960,881 alt≥1 sites, 20× the editors'
count). Parent remains a valid germline mask for it (hom concordance 0.9802); the sample's own
noise floor is what fails.

| sample | study | role | cov | alt≥1 | resid | verdict |
|---|---|---|---:|---:|---:|---|
| nCas9-clone2 | PRJNA1042830 | control | 23.63 | 780,594 | 0.97× | qualified |
| nCas9-clone1 | PRJNA1042830 | control | 23.60 | 854,914 | 1.07× | qualified |
| Parent | PRJNA1042830 | control | 29.29 | 1,135,489 | 1.02× | qualified |
| A3A-Y130F-clone2 | PRJNA1006866 | **editor** | 25.02 | 864,617 | 0.99× | — |
| A3A-Y130F-clone5 | PRJNA1006866 | **editor** | 25.73 | 861,941 | 0.95× | — |
| D10A-clone1 | PRJNA1006866 | control | 31.36 | 6,457,037 | 5.15× | DISQUALIFIED |
| D10A-clone6 | PRJNA1006866 | control | 39.26 | 5,837,065 | 2.96× | DISQUALIFIED |
| P66-background | PRJNA1006866 | control | 37.12 | 16,960,881 | 9.71× | DISQUALIFIED |

**The split is by role, inside one study** — PRJNA1006866 controls 2.96–9.71× (all fail),
PRJNA1006866 editors 0.95–0.99× (both clean), PRJNA1042830 0.97–1.07× (all clean). It is not
"that study is bad."

**Consequence 1.** There is no usable within-study control for the A3A-Y130F arm. The
cross-study nCas9 substitution in §21/§21.2 is the only option, not a convenience.

**Consequence 2, the serious one.** A systematic control-vs-editor split *inside* a single
study points at **processing, not biology**. A 20× gap in alt≥1 at similar depth is not subtle.
If the three controls were prepared or sequenced differently from the two editors, any
editor-vs-control contrast within PRJNA1006866 is confounded by that — and the whole
PRJNA1006866 control set should be treated as unusable rather than merely noisy.

## 23. CORRECTION — gate A0 judged the wrong stratum; the within-study calibrator existed

§21/§21.2/§22 said no usable within-study control exists for the A3A-Y130F arm. **That was
wrong, and it was my gate's fault.** Gate A0 failed calibrators on their **alt≥1** noise floor;
`s7c_editor` analyses **alt≥2**. Over all 23 chromosomes:

| sample | alt≥1 (what the gate judged) | alt≥2 (what the analysis uses) |
|---|---|---|
| D10A-clone1 | 5.13× FAIL | **1.15× PASS** |
| D10A-clone6 | 2.92× FAIL | **0.84× PASS** |
| P66-background | 9.63× FAIL | 4.67× FAIL |

Both D10A clones' entire defect is alt=1 sequencing error — 94.6% / 94.0% singletons against
70–80% in sound samples, mean alt 1.49 / 1.67 against 3.6–4.0 — and the analysis discards that
stratum before computing anything. §22's Consequence 1 is retracted. Background fails both
strata and stays out.

§22's *other* claims survived audit: the extrapolation alternative is excluded
(corr(distance-beyond-fit, residual) = −0.065, non-monotone) and no strand skew exists anywhere
(all eight samples 0.990–1.010).

### 23.1 The within-study arm, run at last — it replicates the null and is better controlled

| comparison | editor MH | calib MH | ed−cal | GC-shift asymmetry |
|---|---:|---:|---:|---|
| clone2 vs D10A-clone6 **within** | 0.777 | 0.961 | **−0.184** | 14.7% vs 10.5% = 4.2pp |
| clone5 vs D10A-clone6 **within** | 0.748 | 0.961 | **−0.213** | 14.1% vs 10.6% = 3.5pp |
| clone2 vs nCas9-clone1 cross | 0.760 | 1.041 | −0.281 | 15.0% vs 1.2% = 13.8pp |
| clone5 vs nCas9-clone1 cross | 0.739 | 1.041 | −0.302 | 14.5% vs 1.2% = 13.3pp |
| **true-null floor** (nCas9-c2 vs c1) | 1.077 | 1.042 | **+0.035** | — |

The GC-shift asymmetry falls from ~14pp to ~4pp — what a same-study calibrator should do — and
the difference shrinks from −0.28/−0.30 to −0.18/−0.21, so **part of the cross-study excess was
study artefact**. Both within-study values remain far on the wrong side of the pre-registered
+0.11 bar and below the measured +0.035 floor. Two clones agreeing to 0.029.

**What does not improve:** burden asymmetry gets *worse* — editor 19,714 specific sites vs
calibrator 171,632 (8.7×, against 4.6× for nCas9). §21.1's caveat stands with more force.
