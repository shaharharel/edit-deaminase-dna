#!/bin/bash
# VA arm driver. VA-clone1 is ~0.3h out and NOTHING would run it -- the haA3A driver was
# written for Y130G only and has already exited. The clone7 driver taught this lesson the
# expensive way: a driver armed without gate A0 fired on a disqualified clone and produced an
# uninterpretable 3-way result. This one gates first, per clone.
#
# ================== PRE-REGISTRATION, WRITTEN BEFORE THE DATA EXISTS ==================
# VA is the OTHER haA3A variant. Y130G, measured earlier today against nCas9-clone1, gave
# GC-adjusted MH ed-cal of +0.318 (clone2) and +0.396 (clone1), against a measured clone-luck
# floor of +0.075 and a pre-registered bar of +0.11.
#
# THE STANDING PRE-REGISTRATION IS UNCHANGED: haA3A was ENGINEERED for near-background
# off-target, so the pre-declared expectation is STILL NO ENRICHMENT. Y130G contradicting it
# does not license me to quietly flip the expectation for VA.
#
# WHAT I AM ADDING IS A DISCRIMINATING PREDICTION, and its value depends entirely on being
# written down now:
#   VA POSITIVE (ed-cal > +0.11, ideally near Y130G's +0.32 to +0.40)
#       -> the effect tracks the haA3A CLASS, not one clone or one library. Two independent
#          variants, four clones total. That is the outcome that would make this worth
#          pursuing seriously.
#   VA NULL (ed-cal inside +/-0.11)
#       -> Y130G's positive is specific to Y130G -- a variant-specific or batch-specific
#          effect. The pre-registration's expectation survives for the class, and the Y130G
#          result becomes a single-variant anomaly needing its own explanation.
#   VA NEGATIVE (ed-cal < -0.11, like A3A-Y130F's -0.18/-0.21)
#       -> the endpoint is producing variant-specific signs with no biological ordering, which
#          is evidence the ENDPOINT is the problem, exactly as the pre-registration warned.
# I commit to reporting whichever of the three lands, and to NOT reinterpreting the categories
# after seeing the number.
# ======================================================================================
BASE=/data/a3a
LOG=$BASE/logs/auto_advance_va.log
exec >> $LOG 2>&1
PY=$HOME/miniconda3/envs/apobec/bin/python
ts(){ date -u +%Y-%m-%dT%H:%M; }

echo "[$(ts)] armed. VA-clone1 and VA-clone2, each processed as it reaches 23/23."
echo "[$(ts)] PRE-REGISTERED EXPECTATION: still NO enrichment (haA3A engineered near-background)."
echo "[$(ts)] DISCRIMINATING PREDICTION recorded above: positive -> class effect; null ->"
echo "[$(ts)]   Y130G-specific; negative -> the endpoint is the problem."
echo "[$(ts)] BAR +0.11 GC-adjusted MH OR. Measured clone-luck floor +0.075. Y130G gave +0.318/+0.396."

done_any=""
for i in $(seq 1 240); do
  for cl in VA-clone1 VA-clone2; do
    case " $done_any " in *" $cl "*) continue;; esac
    n=$(ls $BASE/feat/counts_${cl}_chr*.npz 2>/dev/null | grep -vc partial)
    [ "$n" -lt 23 ] && continue
    echo "[$(ts)] === $cl reached 23/23."
    echo "[$(ts)] GATE A0 on the EDITOR clone itself (the clone7 lesson: gate the editor, not"
    echo "[$(ts)]   just the calibrator -- clone7 was 62x its siblings and nobody checked)."
    if ! A3A_FEAT=$BASE/feat $PY $BASE/qualify_calibrator.py $cl; then
      echo "[$(ts)] *** $cl FAILED GATE A0. NOT running the arm. Recording in BLOCKED.md. ***"
      echo "$(ts) $cl: failed gate A0, VA arm not run for this clone" >> $BASE/logs/BLOCKED.md
      done_any="$done_any $cl"; continue
    fi
    if ! A3A_FEAT=$BASE/feat $PY $BASE/qualify_calibrator.py nCas9-clone1; then
      echo "[$(ts)] *** calibrator nCas9-clone1 failed gate A0. NOT running. ***"
      done_any="$done_any $cl"; continue
    fi
    echo "[$(ts)] === ARM: $cl vs nCas9-clone1 (deaminase-free, same study, Parent-masked)"
    $PY $BASE/s7c_editor.py $cl nCas9-clone1 "" Parent
    echo "[$(ts)] === GC-stratified (read the MH column, not the crude OR)"
    $PY $BASE/qa_editor_gc_B.py $cl nCas9-clone1 || echo "  (gc script signature differs)"
    echo "[$(ts)] === complexity-stratified (check 4's other half)"
    $PY $BASE/qa_complexity_strat.py $cl nCas9-clone1 || true
    echo "[$(ts)] === strand-stratified + read-orientation audit (the bug-4 family)"
    $PY $BASE/qa_strand_hp.py $cl nCas9-clone1 || true
    done_any="$done_any $cl"
  done
  case "$done_any" in *VA-clone1*VA-clone2*|*VA-clone2*VA-clone1*) break;; esac
  sleep 300
done
echo "[$(ts)] done. Read GC-adjusted MH against +0.11, and against Y130G's +0.318/+0.396."
