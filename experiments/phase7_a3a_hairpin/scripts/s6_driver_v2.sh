#!/bin/bash
# Pileup driver v2 — same as v1, plus SELECTIVE BAM RETENTION.
#
# v1 deleted every BAM after its 23 chromosomes, to bound disk. That is right for
# controls, but wrong for the EDITOR arms: if a hairpin result looks striking we
# will want to inspect actual alignments at the top sites (mapping artifacts,
# read orientation, duplicate structure), and regenerating one BAM is ~8h.
#
# Counts keep cov/alt/alt_fwd/alt_rev per site, which covers most follow-ups --
# but not read-level inspection. At 21 GB per BAM against ~525 GB free, retaining
# the primary editor arm is cheap insurance against an expensive re-run.
#
# KEEP_PATTERN samples are retained; everything else is freed as before.
set -uo pipefail
BASE=/mnt/data/a3a
CHROMS=$(seq 1 22; echo X)
PAR=${PAR:-4}
KEEP_PATTERN=${KEEP_PATTERN:-'A3A-Y130F'}
MIN_FREE_KEEP_GB=${MIN_FREE_KEEP_GB:-150}
log(){ echo "[$(date '+%F %T')] $*"; }

log "driver v2 started; retaining BAMs matching '$KEEP_PATTERN' (unless free<${MIN_FREE_KEEP_GB}G)"

while true; do
  for f in $BASE/flags/S2_*_DONE; do
    [ -e "$f" ] || continue
    s=$(basename $f); s=${s#S2_}; s=${s%_DONE}
    bam=$BASE/wgs/${s}.bam
    [ -f "$bam" ] || continue
    [ -f "$BASE/flags/S6_${s}_ALL_DONE" ] && continue
    log "pileup sweep: $s"
    printf '%s\n' $CHROMS | xargs -P $PAR -I{} nice -n 10 \
      $HOME/miniconda3/envs/apobec/bin/python $BASE/s6_pileup.py $s $bam {} \
      >> $BASE/logs/s6_${s}.log 2>&1
    n=$(ls $BASE/feat/counts_${s}_chr*.npz 2>/dev/null | wc -l)
    if [ "$n" -eq 23 ]; then
      touch $BASE/flags/S6_${s}_ALL_DONE
      free=$(df -BG --output=avail /mnt/data | tail -1 | tr -dc '0-9')
      if echo "$s" | grep -qE "$KEEP_PATTERN" && [ "$free" -gt "$MIN_FREE_KEEP_GB" ]; then
        log "$s: 23/23 counts written; BAM RETAINED (editor arm, free ${free}G)"
      else
        rm -f "$bam" "$bam.bai"
        log "$s: 23/23 counts written; BAM freed; disk $(df -h /mnt/data|tail -1|awk '{print $4}')"
      fi
    else
      log "$s: only $n/23 counts -- BAM KEPT for retry"
    fi
  done
  sleep 300
done
