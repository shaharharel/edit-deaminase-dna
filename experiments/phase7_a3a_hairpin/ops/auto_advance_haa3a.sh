#!/bin/bash
# haA3A arm driver — fires when Y130G clone1 AND clone2 both have 23/23 counts.
#
# WHY THIS ARM IS NOW RUNNABLE: SRR26881527-35 (Y130G x2, VA x2, YE1 x2, nCas9 x2, Parent)
# are ALL PRJNA1042830 — one study. The deaminase-free control was never missing, it was
# simply not in node B's queue. nCas9-clone1/2 + Parent counts were transferred from node A
# 2026-08-22 08:35 and verified (0 defects, 236,482,055 sites, strand ~1.01).
#
# STANDING PRE-REGISTRATION, restated before any number exists: haA3A (Y130G/VA) was
# ENGINEERED for near-background off-target. EXPECT NO ENRICHMENT. A null here VALIDATES the
# endpoint. A positive here means the endpoint is broken, not that the editor is active.
# This arm is EXPLORATORY, not confirmatory. Do not restore confirmatory framing afterwards.
#
# BAR (same as node A's, measured, stamped before the run): the clone-luck spread on the
# GC-adjusted MH endpoint is 0.110 and is 5.2x the counting error. An editor effect must
# exceed the calibrator by MORE THAN 0.11 MH OR. Inside 0.11 is null.
BASE=/data/a3a
LOG=$BASE/logs/auto_advance_haa3a.log
exec >> $LOG 2>&1
PY=$HOME/miniconda3/envs/apobec/bin/python
ts(){ date -u +%Y-%m-%dT%H:%M; }

echo "[$(ts)] armed. waiting for Y130G-clone1 and Y130G-clone2 counts (23/23 each)."
echo "[$(ts)] PRE-REGISTERED: expect NO enrichment (haA3A engineered for near-background)."
echo "[$(ts)] BAR: must beat calibrator by >0.11 GC-adjusted MH OR. Inside 0.11 = null."

while :; do
  n1=$(ls $BASE/feat/counts_Y130G-clone1_chr*.npz 2>/dev/null | wc -l)
  n2=$(ls $BASE/feat/counts_Y130G-clone2_chr*.npz 2>/dev/null | wc -l)
  [ "$n1" -ge 23 ] && [ "$n2" -ge 23 ] && break
  sleep 300
done
echo "[$(ts)] both Y130G clones at 23/23. Running the arm."

# 1. within-study VALIDATION first: deaminase-free vs deaminase-free. This must come back
#    near the clone floor. If it does not, nothing downstream is interpretable.
echo "[$(ts)] === VALIDATION: nCas9-clone2 vs nCas9-clone1 (both deaminase-free, same study)"
$PY $BASE/s7c_editor.py nCas9-clone2 nCas9-clone1 '' Parent

# 2. the arm itself, each clone against the real deaminase-free calibrator
for cl in Y130G-clone2 Y130G-clone1; do
  echo "[$(ts)] === ARM: $cl vs nCas9-clone1 (deaminase-free, same study, Parent-masked)"
  $PY $BASE/s7c_editor.py $cl nCas9-clone1 '' Parent
done

# 3. GC-stratified, because every crude OR in this project moves 7-12% under adjustment
echo "[$(ts)] === GC-stratified (crude ORs here are inflated 7-12%; read the MH column)"
for cl in Y130G-clone2 Y130G-clone1; do
  $PY $BASE/qa_editor_gc_B.py $cl nCas9-clone1 || echo "  (gc script signature differs; skipped $cl)"
done

# 4. and the within-family null control that sets this estimator's own floor
echo "[$(ts)] === NULL CONTROL: Y130G-clone1 vs Y130G-clone2 (same editor => clone-luck only)"
$PY $BASE/qa_editor_gc_B.py Y130G-clone1 Y130G-clone2 || true

echo "[$(ts)] done. Read the GC-adjusted MH against the 0.11 bar, not against 1.0."
