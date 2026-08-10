# Cron prompts to recreate on session resume

Both are **session-only** — they die when the Claude session ends and must be
recreated by hand in the new session. Node-side systemd units are unaffected.

---

## MONITOR — every 10 min, cron `4,14,24,34,44,54 * * * *`

```
AUTONOMOUS MONITOR TICK — A3A/PCAWG run. Work on ai-chem (zone us-east1-b,
--tunnel-through-iap). Base dir /mnt/data/a3a (logs/, flags/, pcawg/, wgs/, feat/).
Report only if something changed or needs attention.

1. HEALTH: systemctl is-active for a3a-s2-*, a3a-s6-driver, a3a-reaper.
   df -h /mnt/data (alert if <80G free; /mnt/data/mcgrath_wgs is a RESOLVED null
   and may be deleted for space). Check the instance is not preempted
   (gcloud compute instances describe ai-chem --zone us-east1-b
   --format='value(status)'); if TERMINATED, restart and relaunch every unit whose
   flag file is absent — all scripts are idempotent via /mnt/data/a3a/flags/.
2. FIX FAILURES: tail the relevant log in /mnt/data/a3a/logs/. Diagnose, fix,
   relaunch via sudo systemd-run --unit=<name> --uid=$(id -u) --gid=$(id -g)
   --working-directory=/mnt/data/a3a --setenv=HOME=$HOME. systemctl reset-failed
   before reusing a unit name. Never silently skip a sample — fix it or record why
   in /mnt/data/a3a/logs/BLOCKED.md.
3. KEEP MACHINES BUSY: run pileups for any sample with a DONE flag; advance the
   editor test when an editor arm and its matched calibrator both have counts.
4. IDLE MACHINES: check /proc/<pid>/cmdline for CO-TENANT jobs before stopping
   any node — ai-chem/ai-chem2 are SHARED.
5. Keep /mnt/data/a3a/STATUS.md current.

Discipline: negatives matched on trinucleotide context AND strand; RANDOM baseline
beside every enrichment; deaminase-free clones (nCas9, D10A) plus Parent are the
calibrator for any editor claim, pre-registered bar = editor increment >=3x over
deaminase-free, and the deaminase-free baseline is ~1.4-1.6x NOT 1.0. Never
launder a >1x number into a positive without its random baseline.
```

---

## QA — every 30 min, cron `9,39 * * * *`

```
QA TICK — A3A/PCAWG run. Do this QA YOURSELF, inline. Do NOT spawn subagents —
every subagent in this project has gone idle without returning content, so
delegating means QA silently does not happen.

Node ai-chem: gcloud compute ssh ai-chem --zone us-east1-b --tunnel-through-iap
--command="...". Base /mnt/data/a3a. Python ~/miniconda3/envs/apobec/bin/python.

SIX bugs found so far, all the same family — a strand or coordinate convention
correct locally, then consumed downstream as if universal (see HANDOFF.md §8).
ASSUME A SEVENTH EXISTS.

1. LEAK CHECKS on any new trainset: focal base C for 100% of BOTH classes; ACGT
   distribution at offsets -1 and +1 identical between classes.
2. DERIVE backgrounds, never assume. YTCA background at TCW is 0.6057, not 0.5.
3. RANDOM BASELINE beside every enrichment. Any number equal to 1/base_rate is the
   arithmetic ceiling = leakage, not skill.
4. Signals must survive GC-decile and complexity stratification.
5. Enrichment is the endpoint, not AUROC. Report both; trust enrichment.
6. WATCH FOR SMALL-n ILLUSIONS. The A3A-vs-A3B claim died this way: +0.355 at
   n=53 -> +0.079 at n=97. Report n with every effect; check whether effects
   strengthen or weaken as n grows.
7. EDITOR CLAIMS: no claim without the deaminase-free calibrator (D10A for
   PRJNA1006866, nCas9 for PRJNA1042830, plus Parent/background) and the mandatory
   germline filter. Bar: >=3x increment over deaminase-free, whose baseline is
   ~1.4-1.6x not 1.0. STANDING PRE-REGISTRATION: haA3A (Y130G/VA) was engineered
   for near-background off-target activity — expect NO enrichment. The editor test
   is EXPLORATORY, not confirmatory; do not restore confirmatory framing after
   seeing results.
8. INTEGRITY: /mnt/data/a3a/logs/BLOCKED.md; a3a-reaper active; no S2_*_CLAIMED
   older than 4h without DONE and without a live process.
9. VALIDATE ANY NEW ANALYSIS on a case whose answer is constrained in advance
   (e.g. a deaminase-free clone as pseudo-editor). This is how the 1.4-1.6x
   calibrator baseline was discovered.

If you find a correctness bug, FIX IT and rerun the affected stage immediately.
Append findings to /mnt/data/a3a/STATUS.md. Relay a SHORT summary: what is sound,
what is broken, what you fixed.
```
