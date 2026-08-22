#!/bin/bash
# GC-stratified endpoint B for the clone6 pairs, alongside the main driver.
#
# A SEPARATE driver rather than a fourth edit to auto_advance3.sh, which is armed and waiting.
# Every replacement of a live driver is a chance to break the thing that produces the headline
# result; two read-only analyses on the same trigger cost nothing and risk nothing.
#
# WHY THIS RUNS AT ALL: check 4 had never been applied to the editor arms until 05:35. The
# stem>=6 background spans 2.9x across GC deciles (0.03564 AT-rich -> 0.01248 GC-rich), so any
# shift in the GC composition of the called sites moves the pooled enrichment directly.
# Measured on clone2 vs nCas9: the editor's pooled number carried a 15.0% GC-composition
# effect, the calibrator's 1.2%. An asymmetric confound above the ~12% noise floor is exactly
# what manufactures a spurious between-arm difference.
#
# The banded estimator that is the method of record fixes DEPTH and leaves this untouched.
# Neither substitutes for the other, so both get run.
BASE=/mnt/data/a3a
PY=/home/shaharh_quris_ai/miniconda3/envs/apobec/bin/python
LOG=$BASE/logs/auto_advance_gc.log
exec >> $LOG 2>&1
echo "[$(date -u +%FT%H:%M)] armed; waiting for P66-D10A-clone6 counts (23/23)"
while [ "$(ls $BASE/feat/counts_P66-D10A-clone6_chr*.npz 2>/dev/null | wc -l)" -lt 23 ]; do
  sleep 300
done
echo "[$(date -u +%FT%H:%M)] clone6 complete. GC-stratified endpoint B, within-study pairs:"
for ed in P66-A3A-Y130F-clone2 P66-A3A-Y130F-clone5; do
  echo "=== $ed vs P66-D10A-clone6 (GC-stratified) ==="
  $PY $BASE/qa_editor_gc.py $ed P66-D10A-clone6
done
echo "[$(date -u +%FT%H:%M)] and the deaminase-free pair, as the null control for this estimator:"
echo "=== nCas9-clone2 vs nCas9-clone1 (GC-stratified) ==="
$PY $BASE/qa_editor_gc.py nCas9-clone2 nCas9-clone1
echo "[$(date -u +%FT%H:%M)] done. Read the GC-adjusted MH, not the crude OR."
