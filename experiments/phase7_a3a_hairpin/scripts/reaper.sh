#!/bin/bash
# Clears stale S2_*_CLAIMED flags: a worker killed by preemption between
# touch CLAIMED and touch DONE would otherwise make pop() skip that sample
# forever. A claim with no DONE, older than STALE_MIN and with no live
# aria2c/bwa process for that sample, is released back to the queue.
STALE_MIN=240
F=/mnt/data/a3a/flags
while true; do
  for c in $F/S2_*_CLAIMED; do
    [ -e "$c" ] || continue
    s=$(basename "$c"); s=${s#S2_}; s=${s%_CLAIMED}
    [ -f "$F/S2_${s}_DONE" ] && { rm -f "$c"; continue; }
    if [ -z "$(find "$c" -mmin +$STALE_MIN 2>/dev/null)" ]; then continue; fi
    if pgrep -af 'bwa mem|aria2c' | grep -q "$s"; then continue; fi
    echo "[$(date '+%F %T')] releasing stale claim: $s"
    rm -f "$c"
  done
  sleep 600
done
