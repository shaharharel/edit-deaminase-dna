#!/bin/bash
# Watchdog: retire q1 (old FTP script) the moment its current sample completes.
#
# q1 runs the ORIGINAL s2_queue.sh, which builds ftp:// URLs and has no FAILED
# flag. ENA's FTP now returns "Resource not found" for the PRJNA1006866 paths
# (HTTPS works), and the old script releases the claim on failure and re-claims
# immediately -- it spun 9x in 11s when this first happened. So when
# P66-background finishes, q1 would spin-loop forever on the next sample while
# still reporting "active".
#
# Editing the script in place is unsafe (bash reads it incrementally) and killing
# q1 now would discard ~14h of alignment. So: wait for the DONE flag, then stop q1
# and start a replacement running s2_queue_v2.sh (HTTPS + FAILED guard + backoff).
#
# Idempotent: exits if the replacement already exists.
set -uo pipefail
FLAGS=/mnt/data/a3a/flags
log(){ echo "[$(date '+%F %T')] $*"; }

log "watchdog armed: waiting for S2_P66-background_DONE"
while true; do
  if systemctl is-active --quiet a3a-s2-q5; then
    log "replacement q5 already active; watchdog exiting"; exit 0
  fi
  if [ -f "$FLAGS/S2_P66-background_DONE" ]; then
    log "P66-background DONE -> retiring q1 (old FTP script)"
    sudo systemctl stop a3a-s2-q1 2>/dev/null
    sleep 5
    sudo systemctl reset-failed a3a-s2-q5 2>/dev/null
    sudo systemd-run --unit=a3a-s2-q5 --uid=$(id -u) --gid=$(id -g) \
      --working-directory=/mnt/data/a3a --setenv=HOME="$HOME" \
      --setenv=THREADS=9 --setenv=SORT_THREADS=3 \
      bash -c '/mnt/data/a3a/s2_queue_v2.sh > /mnt/data/a3a/logs/s2_q5.log 2>&1'
    sleep 10
    log "q1=$(systemctl is-active a3a-s2-q1) q5=$(systemctl is-active a3a-s2-q5)"
    exit 0
  fi
  # safety net: if q1 ever starts spinning anyway, catch it by BLOCKED.md growth
  n=$(wc -l < /mnt/data/a3a/logs/BLOCKED.md 2>/dev/null || echo 0)
  if [ "$n" -gt 20 ]; then
    log "BLOCKED.md has $n lines -- q1 may be spinning; retiring it now"
    sudo systemctl stop a3a-s2-q1 2>/dev/null
    sleep 5
    sudo systemctl reset-failed a3a-s2-q5 2>/dev/null
    sudo systemd-run --unit=a3a-s2-q5 --uid=$(id -u) --gid=$(id -g) \
      --working-directory=/mnt/data/a3a --setenv=HOME="$HOME" \
      --setenv=THREADS=9 --setenv=SORT_THREADS=3 \
      bash -c '/mnt/data/a3a/s2_queue_v2.sh > /mnt/data/a3a/logs/s2_q5.log 2>&1'
    exit 0
  fi
  sleep 300
done
