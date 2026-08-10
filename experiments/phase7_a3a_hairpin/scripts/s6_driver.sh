#!/bin/bash
# Runs s6_pileup across all chromosomes for every sample that has finished
# aligning, then frees the BAM once all 23 chromosomes are done.
#
# BAMs are the only real disk consumer (21 GB each, 18 still to come vs ~500 GB
# free). Counts are 13 MB per sample-chromosome, so keeping counts and dropping
# BAMs bounds the footprint. The BAM is only removed after all 23 counts files
# exist -- never on a partial run.
set -uo pipefail
BASE=/mnt/data/a3a
CHROMS=$(seq 1 22; echo X)
PAR=${PAR:-4}
log(){ echo "[$(date '+%F %T')] $*"; }

while true; do
  for f in $BASE/flags/S2_*_DONE; do
    [ -e "$f" ] || continue
    s=$(basename $f); s=${s#S2_}; s=${s%_DONE}
    bam=$BASE/wgs/${s}.bam
    [ -f "$bam" ] || continue
    [ -f "$BASE/flags/S6_${s}_ALL_DONE" ] && continue
    log "pileup sweep: $s"
    printf '%s\n' $CHROMS | xargs -P $PAR -I{} nice -n 10       $HOME/miniconda3/envs/apobec/bin/python $BASE/s6_pileup.py $s $bam {}       >> $BASE/logs/s6_${s}.log 2>&1
    n=$(ls $BASE/feat/counts_${s}_chr*.npz 2>/dev/null | wc -l)
    if [ "$n" -eq 23 ]; then
      touch $BASE/flags/S6_${s}_ALL_DONE
      rm -f "$bam" "$bam.bai"
      log "$s: 23/23 counts written; BAM freed; disk $(df -h /mnt/data|tail -1|awk '{print $4}')"
    else
      log "$s: only $n/23 counts -- BAM KEPT for retry"
    fi
  done
  sleep 300
done
