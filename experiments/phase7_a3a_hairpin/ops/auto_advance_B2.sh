#!/bin/bash
# Fires the pooled 3-clone cross-family analysis the moment Lj-BE clone5 and clone12 have
# counts. They land around 00:30 and 01:15 UTC; waiting for a monitor tick to notice wastes
# time and depends on my being prompted at the right moment.
#
# WHY THREE CLONES AND NOT ONE. The single-clone preview said Lj-BE matches eA3A, which under
# the pre-registration is the "substrate, not enzyme" branch. But the A3A-Y130F arm is the
# standing warning: its apparent stem-8 signal came ENTIRELY from sites shared between two
# clones, and that was invisible until the second clone existed. One clone cannot separate a
# clone-private signal from a shared artefact.
#
# Two things run, in this order:
#   1. the three-way sharing test on the Lj-BE clones (private / 2-of-3 / all-3), which is the
#      instrument that caught the A3A-Y130F artefact
#   2. the banded, depth-matched endpoint B for each Lj-BE clone against pooled eA3A
BASE=/data/a3a
PY=/home/shaharh_quris_ai/miniconda3/envs/apobec/bin/python
LOG=$BASE/logs/auto_advance_B.log
exec >> $LOG 2>&1
echo "[$(date -u +%FT%H:%M)] waiting for Lj-BE clone5 and clone12 counts (23/23 each)"
while [ "$(ls $BASE/feat/counts_P66-Lj-BE-clone5_chr*.npz 2>/dev/null | wc -l)" -lt 23 ] \
   || [ "$(ls $BASE/feat/counts_P66-Lj-BE-clone12_chr*.npz 2>/dev/null | wc -l)" -lt 23 ]; do
  sleep 300
done
echo "[$(date -u +%FT%H:%M)] both landed. Three-way sharing test on the Lj-BE arm:"
$PY $BASE/qa_3way_ljbe.py
echo "[$(date -u +%FT%H:%M)] banded endpoint B, each Lj-BE clone vs eA3A-RL1-clone1:"
for c in P66-Lj-BE-clone5 P66-Lj-BE-clone12; do
  echo "=== $c ==="
  $PY $BASE/s7c_editor.py $c P66-eA3A-RL1-clone1 '' Parent
done
echo "[$(date -u +%FT%H:%M)] done. Results in $LOG"
