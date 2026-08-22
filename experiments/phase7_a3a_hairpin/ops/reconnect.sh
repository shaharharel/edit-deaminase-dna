#!/bin/bash
# ONE-SHOT RECOVERY, to run the moment `gcloud auth login` succeeds.
#
# The backlog accumulated over ~9 blind ticks. Executing it as one idempotent script rather
# than rediscovering it tick by tick removes the risk of dropping one of the four pending
# script pushes, and gets the health sweep in FIRST -- before anything is written -- because
# if a node was preempted, pushing files to it is not the priority.
#
# Idempotent: every write is a copy or an append-once guarded by a marker.
set -uo pipefail
Z="--zone us-east1-b --tunnel-through-iap"
REPO=/Users/shaharharel/Documents/github/edit-deaminase-dna/experiments/phase7_a3a_hairpin
SC=$REPO/scripts
say(){ printf '\n========== %s ==========\n' "$*"; }

say "STEP 1 - HEALTH FIRST (a preempted node changes everything below)"
for n in ai-chem ai-chem2; do
  st=$(gcloud compute instances describe $n $Z --format='value(status)' 2>&1 | tail -1)
  echo "  $n: $st"
  [ "$st" = "TERMINATED" ] && echo "  *** $n TERMINATED -- restart + clear stale claims + relaunch units BEFORE anything else ***"
done

say "STEP 2 - what actually happened while I was blind"
gcloud compute ssh ai-chem $Z --command='
  echo "node time: $(date -u +%FT%TZ)"
  echo "bwa=$(pgrep -c -f "bwa mem") pileup=$(pgrep -cf "[s]6_pileup.py") failed=$(systemctl list-units --state=failed --no-legend|wc -l) free=$(df -h /mnt/data|tail -1|awk "{print \$4}")"
  echo "clone6 counts: $(ls /mnt/data/a3a/feat/counts_P66-D10A-clone6_chr*.npz 2>/dev/null|wc -l)/23"
  echo "--- BLOCKED.md:"; cat /mnt/data/a3a/logs/BLOCKED.md 2>/dev/null || echo "  (empty)"
  echo "--- the critical-path driver (auto_advance3.sh writes auto_advance.log, NOT auto_advance3.log):"
  tail -40 /mnt/data/a3a/logs/auto_advance.log
  echo "--- other node-A drivers:"
  tail -5 /mnt/data/a3a/logs/auto_advance_clone7.log 2>/dev/null
  tail -5 /mnt/data/a3a/logs/auto_advance_gc.log 2>/dev/null
  bash /mnt/data/a3a/progress.sh 2>&1 | tail -6' 2>&1 | tail -70

gcloud compute ssh ai-chem2 $Z --command='
  echo "node time: $(date -u +%FT%TZ)"
  echo "bwa=$(pgrep -c -f "bwa mem") failed=$(systemctl list-units --state=failed --no-legend|wc -l) free=$(df -h /data|tail -1|awk "{print \$4}")"
  echo "Y130G-clone2 counts: $(ls /data/a3a/feat/counts_Y130G-clone2_chr*.npz 2>/dev/null|wc -l)/23"
  echo "--- haA3A arm (gate A0 wired in; PRE-REGISTERED: expect NO enrichment):"
  tail -50 /data/a3a/logs/auto_advance_haa3a.log 2>/dev/null
  bash /data/a3a/progress.sh 2>&1 | tail -4' 2>&1 | tail -60

say "STEP 3 - push the four generator fixes made while blind"
gcloud compute scp $SC/qa_donor_jackknife.py ai-chem2:/data/a3a/qa_donor_jackknife.py $Z 2>&1|tail -1
for f in pcawg_stage2.py pcawg_stage2b.py pcawg_build.py; do
  gcloud compute scp $SC/$f ai-chem2:/data/a3a/$f $Z 2>&1|tail -1
done
gcloud compute ssh ai-chem2 $Z --command='
  cd /data/a3a
  for f in qa_donor_jackknife.py pcawg_stage2.py pcawg_stage2b.py pcawg_build.py; do
    ~/miniconda3/envs/apobec/bin/python -c "compile(open(\"$f\").read(),\"$f\",\"exec\"); print(\"  $f compiles\")" 2>&1|tail -1
  done
  echo "  guards present:"
  grep -lc "A3A_ALLOW_V1_RANKING\|A3A_ALLOW_BROKEN_TETRA" pcawg_stage2.py pcawg_stage2b.py pcawg_build.py 2>/dev/null
  grep -c "np.ptp" qa_donor_jackknife.py'

say "STEP 4 - append the pending STATUS entries (append-once)"
for f in status56.md status57.md; do
  [ -e "$SCRATCH/$f" ] || continue
done
echo "  (scp + append status56.md and status57.md from the scratchpad, guarded by a grep for"
echo "   their headline strings so a rerun does not duplicate them)"

say "STEP 5 - THEN resume the normal tick"
echo "  read the driver output above against the PRE-REGISTERED bars:"
echo "    node A  A3A-Y130F vs D10A-clone6 must beat the calibrator by >0.11 GC-adjusted MH OR"
echo "    node B  haA3A expects NO enrichment; a null VALIDATES the endpoint"
echo "  and check idle machines against /proc/<pid>/cmdline for CO-TENANTS before any stop."
