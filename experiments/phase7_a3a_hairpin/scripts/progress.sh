#!/bin/bash
# TRUE s2 progress. Replaces the read-count metric, which is WRONG.
#
# WHY: logs/s2_<sample>_bwa.log is APPEND-ONLY and was created 2026-08-10. Each log
# therefore holds an aborted 10-day-old run PLUS the current run, concatenated with no
# banner between them. Summing 'read N sequences' over the whole file adds the dead
# run's reads to the live one and roughly DOUBLES the apparent progress -- which is how
# clone2 reported "103%" while actually sitting near 50%.
# Tell: two distinct batch sizes in one log. 600000 = the dead 2026-08-10 generation,
# 466668 = the current -t 7 generation. nCas9-clone2, which completed under the OLD
# generation, shows only 600000 -- so the batch size is a reliable generation marker.
#
# The authoritative signal is the byte offset of bwa's own open FASTQ fd: it carries no
# history, and the two mate files agree to ~0.1% which is a built-in consistency check.
CUR_BATCH=466668
printf '%-26s %7s %7s %12s %12s\n' SAMPLE MATE1 MATE2 CUR_READS ETA
for p in $(pgrep -f 'bwa mem'); do
  s=$(tr '\0' '\n' < /proc/$p/cmdline | grep -o 'ID:[^\\]*' | head -1 | cut -d: -f2)
  [ -z "$s" ] && continue
  pcts=();
  for fd in /proc/$p/fd/*; do
    t=$(readlink "$fd" 2>/dev/null); case "$t" in *fastq.gz|*fq.gz) ;; *) continue;; esac
    n=$(basename "$fd"); pos=$(awk '/^pos:/{print $2}' /proc/$p/fdinfo/$n 2>/dev/null)
    sz=$(stat -c%s "$t" 2>/dev/null)
    [ -n "$pos" ] && [ -n "$sz" ] && [ "$sz" -gt 0 ] && pcts+=("$(awk -v a="$pos" -v b="$sz" 'BEGIN{printf "%.2f", a*100/b}')")
  done
  L=/mnt/data/a3a/logs/s2_${s}_bwa.log
  cur=$(grep -ao 'read [0-9]* sequences' "$L" 2>/dev/null | awk -v c=$CUR_BATCH '$2==c{n++} END{print (n+0)*c}')
  el=$(ps -o etimes= -p $p | tr -d ' ')
  pct=${pcts[0]:-0}
  eta=$(awk -v e="$el" -v p="$pct" 'BEGIN{ if(p>0.5) printf "%.1fh", (e*(100-p)/p)/3600; else print "?" }')
  printf '%-26s %6s%% %6s%% %12d %12s\n' "$s" "${pcts[0]:-?}" "${pcts[1]:-?}" "$cur" "$eta"
done
