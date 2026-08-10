#!/bin/bash
# Queue-driven WGS worker. Pops SAMPLE:SRR lines from a shared queue under flock,
# so N workers stay busy and new datasets can be appended to the queue at any time
# without restarting anything. Idempotent via flags/.
#
# Usage: THREADS=9 s2_queue.sh
set -uo pipefail

BASE=/mnt/data/a3a
WGS=$BASE/wgs
FLAGDIR=$BASE/flags
LOGS=$BASE/logs
QUEUE=$BASE/queue.txt
QLOCK=$BASE/queue.lock
REF=/mnt/data/ref/hg19.fa
THREADS=${THREADS:-9}
SORT_THREADS=${SORT_THREADS:-3}
SORT_MEM=${SORT_MEM:-2G}
MIN_FREE_GB=${MIN_FREE_GB:-120}

export PATH="$HOME/miniconda3/envs/bio/bin:$PATH"
mkdir -p "$WGS" "$FLAGDIR" "$LOGS"
log() { echo "[$(date '+%F %T')] $*"; }

pop() {   # atomically take the first not-yet-done line off the queue
  flock "$QLOCK" bash -c '
    Q='"$QUEUE"'; F='"$FLAGDIR"'
    [ -s "$Q" ] || exit 1
    while read -r line; do
      [ -z "$line" ] && continue
      s=${line%%:*}
      if [ ! -f "$F/S2_${s}_DONE" ] && [ ! -f "$F/S2_${s}_CLAIMED" ] && [ ! -f "$F/S2_${s}_FAILED" ]; then
        touch "$F/S2_${s}_CLAIMED"; echo "$line"; exit 0
      fi
    done < "$Q"
    exit 1'
}

ena_url() {
  curl -s --max-time 90 --retry 3 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" \
    | tr '\t' '\n' | grep 'ftp.sra.ebi.ac.uk' | tr ';' '\n' | sed 's#^#https://#'
}

while true; do
  spec=$(pop) || { log "queue empty, worker exiting"; break; }
  SAMPLE=${spec%%:*}; SRR=${spec##*:}
  BAM=$WGS/${SAMPLE}.bam

  free_gb=$(df -BG --output=avail /mnt/data | tail -1 | tr -dc '0-9')
  if [ "$free_gb" -lt "$MIN_FREE_GB" ]; then
    log "LOW DISK (${free_gb}G < ${MIN_FREE_GB}G) — releasing $SAMPLE and pausing 15m"
    rm -f "$FLAGDIR/S2_${SAMPLE}_CLAIMED"; sleep 900; continue
  fi

  log "=== $SAMPLE ($SRR) start (free ${free_gb}G)"
  fail=0
  for i in 1 2; do
    F=$WGS/${SRR}_${i}.fastq.gz
    if [ -f "$F" ] && gzip -t "$F" 2>/dev/null; then log "  read$i cached"; continue; fi
    URL=$(ena_url "$SRR" | sed -n "${i}p")
    if [ -z "$URL" ]; then log "  FAIL: no ENA URL for $SRR read$i"; fail=1; break; fi
    aria2c -x16 -s16 -c --file-allocation=none --summary-interval=300 \
           -d "$WGS" -o "${SRR}_${i}.fastq.gz" "$URL" >> "$LOGS/s2_${SAMPLE}_dl.log" 2>&1
    gzip -t "$F" 2>/dev/null || { log "  FAIL: $F corrupt"; fail=1; break; }
  done
  if [ $fail -ne 0 ]; then
    echo "$SAMPLE ($SRR): download failed $(date)" >> "$LOGS/BLOCKED.md"
    touch "$FLAGDIR/S2_${SAMPLE}_FAILED"
    rm -f "$FLAGDIR/S2_${SAMPLE}_CLAIMED"; sleep 60; continue
  fi

  log "  bwa mem -> sorted BAM (${THREADS}t)"
  bwa mem -t $THREADS -M \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
      "$REF" "$WGS/${SRR}_1.fastq.gz" "$WGS/${SRR}_2.fastq.gz" \
      2>> "$LOGS/s2_${SAMPLE}_bwa.log" \
    | samtools sort -@ $SORT_THREADS -m $SORT_MEM -T "$WGS/tmp_${SAMPLE}" \
      -o "$BAM" - 2>> "$LOGS/s2_${SAMPLE}_bwa.log"

  if [ ! -s "$BAM" ] || ! samtools quickcheck "$BAM"; then
    log "  FAIL: BAM bad for $SAMPLE"
    echo "$SAMPLE: alignment failed $(date)" >> "$LOGS/BLOCKED.md"
    rm -f "$BAM" "$FLAGDIR/S2_${SAMPLE}_CLAIMED"; continue
  fi
  samtools index -@ 8 "$BAM"
  samtools idxstats "$BAM" > "$WGS/${SAMPLE}.idxstats"
  samtools flagstat -@ 8 "$BAM" > "$WGS/${SAMPLE}.flagstat"
  rm -f "$WGS/${SRR}_1.fastq.gz" "$WGS/${SRR}_2.fastq.gz"
  touch "$FLAGDIR/S2_${SAMPLE}_DONE"
  log "=== $SAMPLE DONE ($(du -h "$BAM"|cut -f1)); free $(df -h /mnt/data|tail -1|awk '{print $4}')"
done
