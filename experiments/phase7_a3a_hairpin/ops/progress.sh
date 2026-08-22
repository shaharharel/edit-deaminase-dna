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
# CUR_BATCH is DERIVED per sample below, not hardcoded: it is 466668 at 7 threads and
# 933334 at 14, so a fixed value silently zeroes the column on the other node.
printf '%-26s %6s %6s %12s %11s %9s %9s\n' SAMPLE MATE1 MATE2 CUR_READS IMPLIED ETA_AVG LAST_OUT
for p in $(pgrep -f 'bwa mem'); do
STATE=${STATE:-/tmp/progress_state.txt}   # rolling offsets for the RECENT-rate column
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
  # last matching line = the LIVE generation; these logs are append-only and can hold a
  # preempted run's batches too, so most-frequent would be wrong once a dead run ran longer
  cb=$(grep -ao 'read [0-9]* sequences' "$L" 2>/dev/null | awk '{b=$2} END{print b+0}')
  cur=$(grep -ao 'read [0-9]* sequences' "$L" 2>/dev/null | awk -v c="$cb" '$2==c{n++} END{print (n+0)*c}')
  el=$(ps -o etimes= -p $p | tr -d ' ')
  pct=${pcts[0]:-0}
  eta=$(awk -v e="$el" -v p="$pct" 'BEGIN{ if(p>0.5) printf "%.1fh", (e*(100-p)/p)/3600; else print "?" }')
  # implied total = counter extrapolated by the file fraction. If a sample was restarted
  # with the same thread count, the log holds two indistinguishable generations and the
  # counter double-counts -- which shows up here as an implausible implied total.
  imp=$(awk -v c="$cur" -v p="$pct" 'BEGIN{ if(p>0.5) printf "%.0f", c*100/p; else print 0 }')
  # RECENT rate, from the previous invocation's stamp. The lifetime ETA above cannot see a
  # slowdown; this can. STATE holds "sample epoch pct" lines.
  now=$(date +%s)
  # QUIET TIME, not a rate. Batches take 8-25 min and this is sampled every ~10 min, so a
  # two-point rate is dominated by batch-boundary placement: it produced 30.1h, 42.7h, 30.0h
  # and 103.9h tonight, all false, and zero true positives. Quiet time is a direct measurement
  # and is what every intervention decision actually used.
  reta=$(awk -v a="$(( now - $(stat -c %Y "$L" 2>/dev/null || echo $now) ))" \
         'BEGIN{ if(a<90) printf "active"; else printf "%.0fm quiet", a/60 }')
  echo "$s $now $pct" >> "$STATE"
  printf '%-26s %6s%% %6s%% %12d %10sM %9s %9s\n' "$s" "${pcts[0]:-?}" "${pcts[1]:-?}" "$cur" \
    "$(awk -v i="$imp" 'BEGIN{printf "%.0f", i/1e6}')" "$eta" "$reta"
done

# --- samples past alignment and still in sort/merge. bwa has exited so they have no FASTQ
# descriptor and would otherwise be invisible here, which reads as "finished" or "gone".
for pid in $(pgrep -f "[s]amtools sort" 2>/dev/null); do
  smp=$(tr "\0" " " < /proc/$pid/cmdline 2>/dev/null | grep -o "tmp_[A-Za-z0-9-]*" | head -1)
  smp=${smp#tmp_}
  [ -z "$smp" ] && continue
  # a live samtools sort runs THROUGHOUT alignment, spilling as it goes. Only when bwa for
  # this sample has exited is the sort actually in its final merge.
  pgrep -f "bwa mem.*SM:${smp}\b" > /dev/null 2>&1 && continue
  # the sort process is as old as the ALIGNMENT -- it streams bwa's output the whole way.
  # The merge began when bwa exited, i.e. at the bwa log's last write.
  B=${BASE:-/mnt/data/a3a}
  lg=$B/logs/s2_${smp}_bwa.log
  mstart=$(stat -c %Y "$lg" 2>/dev/null || echo 0)
  mmin=$(( ( $(date +%s) - mstart ) / 60 ))
  spills=$(ls $B/wgs/tmp_${smp}.*.bam 2>/dev/null | wc -l)
  gb=$(awk -v b="$(stat -c %s $B/wgs/${smp}.bam 2>/dev/null || echo 0)" 'BEGIN{printf "%.1fG out", b/1e9}')
  printf '%-26s %6s %6s %12s %11s %9s %9s\n' "$smp" "MERGING" "-" "$spills spills" "$gb" \
    "${mmin}m merging" "-"
done
