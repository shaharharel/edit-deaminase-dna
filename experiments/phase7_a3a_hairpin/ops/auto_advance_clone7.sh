#!/bin/bash
# Fires the THREE-WAY test on the A3A-Y130F arm when clone7's counts complete (~6 h).
#
# A separate driver rather than an edit to auto_advance3.sh, which is running. Editing a live
# bash script re-reads it by byte offset; I did exactly that an hour ago by habit and had to
# kill and relaunch the driver to repair it. New file, every time.
#
# WHY THIS TEST. The arm was downgraded on a TWO-clone recurrence result: 3,000 sites shared
# against 1.8 expected (1642x). Two clones cannot separate "systematic artefact" from "shared
# sub-clonal ancestry" -- both produce 2-way sharing. Three can: under independence the counts
# fall off as A*B/E for pairs and A*B*C/E^2 for triples, and the observed ratio between those
# distinguishes the two.
#
# PRE-REGISTERED, before clone7 exists:
#   - artefact  -> the ALL-3 class is large AND hairpin-rich relative to private
#   - ancestry  -> ALL-3 large but NOT hairpin-rich (what eA3A showed: private 1.371x,
#                  all-three 1.110x)
#   - neither   -> 3-way sharing near its chance expectation
# And the standing caveat from tonight: read the RECURRENCE first. The shared-vs-private
# hairpin split was underpowered at n_hp=8 in the 2-clone version (Fisher p=0.520).
BASE=/mnt/data/a3a
PY=/home/shaharh_quris_ai/miniconda3/envs/apobec/bin/python
LOG=$BASE/logs/auto_advance_clone7.log
exec >> $LOG 2>&1
echo "[$(date -u +%FT%H:%M)] armed; waiting for P66-A3A-Y130F-clone7 counts (23/23)"
while [ "$(ls $BASE/feat/counts_P66-A3A-Y130F-clone7_chr*.npz 2>/dev/null | wc -l)" -lt 23 ]; do
  sleep 300
done
echo "[$(date -u +%FT%H:%M)] clone7 complete. Noise-floor qualification first (within-study):"
$PY $BASE/qualify_calibrator.py P66-A3A-Y130F-clone7
echo "[$(date -u +%FT%H:%M)] three-way sharing test on the A3A-Y130F arm:"
$PY $BASE/qa_3way_a3a.py
echo "[$(date -u +%FT%H:%M)] and the clonal-power check on the enlarged arm:"
$PY $BASE/qa_clonal_power.py
echo "[$(date -u +%FT%H:%M)] done."
