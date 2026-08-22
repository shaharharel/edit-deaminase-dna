#!/bin/bash
# Fires the editor analysis when the clean calibrator's counts complete.
#
# v2 ADDS THE RECURRENCE-FILTER COMPARISON, which has been a recorded pending task since the
# filter was written: "re-run with and without A3A_EXCLUDE_RECURRENT so the filter's effect is
# MEASURED, not assumed." The filter exists because two independent A3A-Y130F clones share
# 3,000 editor-specific sites against 1.8 expected by chance, and splitting on that sharing
# showed the marginal stem-8 signal came entirely from the shared fraction (1.326x shared vs
# 1.026x private). Running only the filtered version would bake that judgement in; running
# only the unfiltered version ignores it. Both, side by side, is the honest form.
#
# Order is deliberate and unchanged: the cross-study gate FIRST -- if Parent has drifted from
# the P66 lineage it must not be the germline mask and everything downstream is wrong -- then
# the deaminase-free validation, then the editor arms.
#
# D10A-clone1 is deliberately absent: 6.46M alt>=1 sites with 94.8% below VAF 0.05
# disqualified it. clone6 replaces it.
BASE=/mnt/data/a3a
PY=/home/shaharh_quris_ai/miniconda3/envs/apobec/bin/python
LOG=$BASE/logs/auto_advance.log
exec >> $LOG 2>&1
echo "[$(date -u +%FT%H:%M)] v2 armed; waiting for P66-D10A-clone6 counts (23/23)"
while [ "$(ls $BASE/feat/counts_P66-D10A-clone6_chr*.npz 2>/dev/null | wc -l)" -lt 23 ]; do
  sleep 300
done
echo "[$(date -u +%FT%H:%M)] clone6 complete. Step A: cross-study gate."
$PY $BASE/s8_xstudy.py Parent P66-D10A-clone6
echo "[$(date -u +%FT%H:%M)] Step B: deaminase-free validation (expect ~1.0 per band)."
$PY $BASE/s7c_editor.py nCas9-clone2 nCas9-clone1 '' Parent
echo "[$(date -u +%FT%H:%M)] Step C: editor arms vs the CLEAN calibrator, UNFILTERED."
for ed in P66-A3A-Y130F-clone2 P66-A3A-Y130F-clone5; do
  echo "=== $ed vs P66-D10A-clone6  [recurrence filter OFF] ==="
  $PY $BASE/s7c_editor.py $ed P66-D10A-clone6 '' Parent
done
echo "[$(date -u +%FT%H:%M)] Step D: the same two arms with the cross-clone recurrence filter ON."
echo "=== P66-A3A-Y130F-clone2, excluding sites also specific in clone5 ==="
A3A_EXCLUDE_RECURRENT=P66-A3A-Y130F-clone5 $PY $BASE/s7c_editor.py \
  P66-A3A-Y130F-clone2 P66-D10A-clone6 '' Parent
echo "=== P66-A3A-Y130F-clone5, excluding sites also specific in clone2 ==="
A3A_EXCLUDE_RECURRENT=P66-A3A-Y130F-clone2 $PY $BASE/s7c_editor.py \
  P66-A3A-Y130F-clone5 P66-D10A-clone6 '' Parent
echo "[$(date -u +%FT%H:%M)] done. Compare Step C against Step D: if the banded MH values"
echo "  move materially when the shared sites are dropped, the arm's signal lives in them."
