#!/bin/bash
# S2 — PRJNA1042830 clonal WGS: download -> align -> sort -> index.
# "Engineered APOBEC3A deaminase for highly accurate cytosine base editing"
# Arms: Y130G (A3A-class editor), VA, YE1 (rAPOBEC1), nCas9 (deaminase-free), Parent.
#
# Idempotent: per-sample flags in $FLAGDIR. Safe to re-run after preemption.
# Disk-frugal: FASTQ deleted once the BAM is verified.
#
# Usage: s2_wgs.sh SAMPLE:SRR [SAMPLE:SRR ...]

set -uo pipefail

BASE=/mnt/data/a3a
WGS=$BASE/wgs
FLAGDIR=$BASE/flags
LOGS=$BASE/logs
REF=/mnt/data/ref/hg19.fa
THREADS=${THREADS:-9}        # 3 concurrent workers x 9 threads beats 1 x 28:
SORT_THREADS=${SORT_THREADS:-3}   # bwa scales sublinearly past ~16 threads
SORT_MEM=${SORT_MEM:-2G}

export PATH="$HOME/miniconda3/envs/bio/bin:$PATH"
mkdir -p "$WGS" "$FLAGDIR" "$LOGS"

log() { echo "[$(date '+%F %T')] $*"; }

ena_url() {  # $1 = SRR ; echoes the two fastq URLs (asked, not constructed --
             # ENA's subdirectory rule is not simply the last 3 digits)
  local srr=$1
  curl -s --max-time 90 --retry 3 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=fastq_ftp&format=tsv" \
    | tr '\t' '\n' | grep 'ftp.sra.ebi.ac.uk' | tr ';' '\n' | sed 's#^#ftp://#'
}

for spec in "$@"; do
  SAMPLE=${spec%%:*}
  SRR=${spec##*:}
  BAM=$WGS/${SAMPLE}.bam

  if [ -f "$FLAGDIR/S2_${SAMPLE}_DONE" ]; then
    log "$SAMPLE already aligned, skipping"
    continue
  fi

  log "=== $SAMPLE ($SRR) starting"

  # ---- download (aria2c resumes; verify with gzip -t, NOT file size:
  #      aria2c fallocates so stat/du show full size immediately)
  for i in 1 2; do
    F=$WGS/${SRR}_${i}.fastq.gz
    if [ -f "$F" ] && gzip -t "$F" 2>/dev/null; then
      log "  $SRR read$i already complete"
      continue
    fi
    URL=$(ena_url "$SRR" | sed -n "${i}p")
    log "  downloading $URL"
    aria2c -x16 -s16 -c --file-allocation=none --summary-interval=120 \
           -d "$WGS" -o "${SRR}_${i}.fastq.gz" "$URL" \
           >> "$LOGS/s2_${SAMPLE}_dl.log" 2>&1
    if ! gzip -t "$F" 2>/dev/null; then
      log "  FAIL: $F corrupt after download, skipping $SAMPLE"
      continue 2
    fi
  done

  # ---- align
  log "  bwa mem -> sorted BAM"
  bwa mem -t $THREADS -M \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
      "$REF" "$WGS/${SRR}_1.fastq.gz" "$WGS/${SRR}_2.fastq.gz" \
      2>> "$LOGS/s2_${SAMPLE}_bwa.log" \
    | samtools sort -@ $SORT_THREADS -m $SORT_MEM -T "$WGS/tmp_${SAMPLE}" \
      -o "$BAM" - 2>> "$LOGS/s2_${SAMPLE}_bwa.log"

  if [ ! -s "$BAM" ]; then
    log "  FAIL: empty BAM for $SAMPLE"
    continue
  fi

  samtools index -@ 8 "$BAM"
  if ! samtools quickcheck "$BAM"; then
    log "  FAIL: quickcheck failed for $SAMPLE"
    continue
  fi

  # ---- record depth so downstream can gate on coverage (the artifact that
  #      dominated every prior null in this project)
  samtools idxstats "$BAM" > "$WGS/${SAMPLE}.idxstats"
  samtools flagstat -@ 8 "$BAM" > "$WGS/${SAMPLE}.flagstat"

  rm -f "$WGS/${SRR}_1.fastq.gz" "$WGS/${SRR}_2.fastq.gz"
  touch "$FLAGDIR/S2_${SAMPLE}_DONE"
  log "=== $SAMPLE DONE ($(du -h "$BAM" | cut -f1)); free: $(df -h /mnt/data | tail -1 | awk '{print $4}')"
done

log "S2 batch complete"
