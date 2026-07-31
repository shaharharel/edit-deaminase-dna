#!/usr/bin/env bash
# Path A v9 — cross-deaminase multi-clone WGS acquisition for the DNA editing index.
# Study: PRJNA553240 (Doman 2020). Cohort: BE4 (rAPOBEC1) x8, YE1-BE4 (engineered) x8,
#        nCas9 (no-deaminase control) x7, Parent (untreated baseline) x1.
#
# Per clone: download paired FASTQ -> bwa mem | samtools sort (RAM-heavy, fixes Path A I/O wall)
#            -> GATK MarkDuplicates -> samtools mpileup at CDS C positions (UNCAPPED depth,
#            Levanon editing-index principle) -> keep only the slim .mpileup, delete FASTQ+BAM.
#
# Idempotent: skips any clone whose final .mpileup already exists (>100 MB). Safe to re-run.
# Sequential by design: one clone at a time avoids the disk I/O contention that killed v8.
set -uo pipefail

DATA=/mnt/data
BIN=$DATA/tools/miniconda/envs/doman/bin
export PATH="$BIN:$DATA/tools/miniconda/bin:$PATH"   # gatk wrapper needs python/java on PATH
SAMTOOLS=$BIN/samtools
BWA=$BIN/bwa
GATK=$BIN/gatk
REF=$DATA/ref/hg19.fa
BED=$DATA/dna_features/index/positions_hg19.bed
GOLD=$DATA/gold_index
LOGDIR=$DATA/logs
LOG=$LOGDIR/pathA_v9.log
LOCK=$DATA/pathA_v9.lock

THREADS="${THREADS:-24}"        # override per machine (32-core->30, 16-core->16)
SORT_MEM="${SORT_MEM:-3G}"      # THREADS x SORT_MEM must stay well under RAM
GCS_OUT="${GCS_OUT:-}"          # if set (gs://.../mpileups), each finished mpileup is uploaded there (durable)
MPILEUP_D=100000       # effectively uncapped for ~150x CDS coverage (no hard coverage cap)
MIN_FREE_GB=350        # worst-case transient peak: fastq + sort tmp + sorted.bam + dedup.bam for a ~140x clone
SORTTMP=$DATA/tmp      # on sdb (813G): capacity-safe; -m 3G keeps writes large/sequential so throughput holds

mkdir -p "$LOGDIR" "$SORTTMP"
exec 9>"$LOCK"
if ! flock -n 9; then echo "another pathA_v9 instance is running; exiting" >&2; exit 0; fi

log(){ echo "[$(date +%Y-%m-%d_%H:%M:%S)] $*" | tee -a "$LOG"; }

# label  ->  ENA fastq prefix (append _1.fastq.gz / _2.fastq.gz). Verified from ENA filereport.
declare -A URL=(
  [BE4_clone1]=SRR104/004/SRR10413104/SRR10413104
  [BE4_clone2]=SRR104/003/SRR10413103/SRR10413103
  [BE4_clone3]=SRR104/092/SRR10413092/SRR10413092
  [BE4_clone4]=SRR104/087/SRR10413087/SRR10413087
  [BE4_clone5]=SRR104/086/SRR10413086/SRR10413086
  [BE4_clone6]=SRR104/085/SRR10413085/SRR10413085
  [BE4_clone7]=SRR104/084/SRR10413084/SRR10413084
  [BE4_clone8]=SRR104/083/SRR10413083/SRR10413083
  [Parent_WGS]=SRR104/088/SRR10413088/SRR10413088
  [YE1-BE4_clone1]=SRR104/082/SRR10413082/SRR10413082
  [YE1-BE4_clone2]=SRR104/081/SRR10413081/SRR10413081
  [YE1-BE4_clone3]=SRR104/002/SRR10413102/SRR10413102
  [YE1-BE4_clone4]=SRR104/001/SRR10413101/SRR10413101
  [YE1-BE4_clone5]=SRR104/000/SRR10413100/SRR10413100
  [YE1-BE4_clone6]=SRR104/099/SRR10413099/SRR10413099
  [YE1-BE4_clone7]=SRR104/098/SRR10413098/SRR10413098
  [YE1-BE4_clone8]=SRR104/097/SRR10413097/SRR10413097
  [nCas9_clone1]=SRR104/096/SRR10413096/SRR10413096
  [nCas9_clone2]=SRR104/095/SRR10413095/SRR10413095
  [nCas9_clone3]=SRR104/094/SRR10413094/SRR10413094
  [nCas9_clone4]=SRR104/093/SRR10413093/SRR10413093
  [nCas9_clone5]=SRR104/091/SRR10413091/SRR10413091
  [nCas9_clone6]=SRR104/090/SRR10413090/SRR10413090
  [nCas9_clone7]=SRR104/089/SRR10413089/SRR10413089
)

# Processing order: reuse the two existing dedup BAMs first (cheap mpileup-only, makes the
# index settings consistent), then interleave BE4/YE1/nCas9 so a usable multi-clone set
# (>=4 BE4, >=3 YE1, >=3 nCas9) lands as early as possible.
ORDER=(
  BE4_clone1 Parent_WGS
  BE4_clone2 YE1-BE4_clone1 nCas9_clone1
  BE4_clone3 YE1-BE4_clone2 nCas9_clone2
  BE4_clone4 YE1-BE4_clone3 nCas9_clone3
  BE4_clone5 YE1-BE4_clone4 nCas9_clone4
  BE4_clone6 YE1-BE4_clone5 nCas9_clone5
  BE4_clone7 YE1-BE4_clone6 nCas9_clone6
  BE4_clone8 YE1-BE4_clone7 nCas9_clone7
  YE1-BE4_clone8
)
# Allow restricting to a subset via args (e.g. for splitting across machines): pathA_v9.sh BE4_clone5 ...
if [ "$#" -gt 0 ]; then ORDER=("$@"); fi

free_gb(){ df -BG --output=avail "$DATA" | tail -1 | tr -dc '0-9'; }

process_clone(){
  local L=$1
  local FINAL=$GOLD/${L}_samtools.mpileup
  if [ -s "$FINAL" ] && [ -f "${FINAL}.done" ]; then
    log "[skip] $L — verified mpileup present ($(du -h "$FINAL"|cut -f1), $(cat "${FINAL}.done") lines)"; return 0
  fi
  rm -f "$FINAL" "${FINAL}.tmp" "${FINAL}.done"   # any leftover from an interrupted run is unverified
  local BAMDIR=$DATA/bam/$L
  local DEDUP=$BAMDIR/${L}.dedup.bam
  local DOWNLOADED=0

  if [ -s "$DEDUP" ]; then
    log "[$L] reusing existing dedup BAM ($(du -h "$DEDUP"|cut -f1)) — mpileup only"
  else
    if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
      log "[$L] ABORT-CLONE: only $(free_gb)G free on $DATA (<${MIN_FREE_GB}G)"; return 1
    fi
    local pre="https://ftp.sra.ebi.ac.uk/vol1/fastq/${URL[$L]}"
    local SRR="${URL[$L]##*/}"
    local FQDIR=$DATA/fastq/$L
    rm -rf "$FQDIR" "$BAMDIR"; mkdir -p "$FQDIR" "$BAMDIR"
    DOWNLOADED=1
    # Best-effort expected sizes / read count from ENA, for integrity checks below
    local META EXP_BYTES READ_COUNT
    # ENA always prepends run_accession as col1, so fastq_bytes=$2 (a;b list), read_count=$3.
    # Use printf %d (not print s+0) so big sums never come back in scientific notation.
    META=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_bytes,read_count&format=tsv" 2>>"$LOG" | tail -1)
    # %.0f (not %d): mawk's %d is 32-bit and overflows on ~64 GB byte sums; doubles are exact here.
    EXP_BYTES=$(echo "$META" | awk -F'\t' '{n=split($2,a,";"); s=0; for(i=1;i<=n;i++) s+=a[i]; printf "%.0f", s}')
    READ_COUNT=$(echo "$META" | awk -F'\t' '{printf "%.0f", $3}')
    log "[$L] $SRR expected_bytes=${EXP_BYTES} read_count=${READ_COUNT}"
    log "[$L] downloading ${URL[$L]}_{1,2}.fastq.gz"
    local i
    for i in 1 2; do
      # aria2c multi-connection (~20x single-stream curl; ENA allows parallel ranges). Resumes via -c.
      aria2c -x16 -s16 -j1 -c --max-tries=20 --retry-wait=15 --connect-timeout=60 --timeout=120 \
        --auto-file-renaming=false --allow-overwrite=true --console-log-level=warn --summary-interval=0 \
        -d "$FQDIR" -o "${L}_${i}.fastq.gz" "${pre}_${i}.fastq.gz" 2>>"$LOG" \
        || { log "[$L] FAIL download R$i"; rm -rf "$FQDIR" "$BAMDIR"; return 1; }
    done
    local GOT_BYTES; GOT_BYTES=$(stat -c %s "$FQDIR/${L}_1.fastq.gz" "$FQDIR/${L}_2.fastq.gz" | awk '{s+=$1} END{print s}')
    if [ "${EXP_BYTES:-0}" -gt 0 ] 2>/dev/null; then
      if [ "$GOT_BYTES" -lt $(( EXP_BYTES * 98 / 100 )) ]; then
        log "[$L] FAIL download incomplete: got=$GOT_BYTES expected=$EXP_BYTES"; rm -rf "$FQDIR" "$BAMDIR"; return 1
      fi
    elif [ "$GOT_BYTES" -lt 2000000000 ]; then
      log "[$L] FAIL download too small: got=$GOT_BYTES (no ENA size to compare)"; rm -rf "$FQDIR" "$BAMDIR"; return 1
    fi
    log "[$L] downloaded $(du -sh "$FQDIR"|cut -f1); bwa mem | sort (-@ $THREADS -m $SORT_MEM)"
    local ST="$SORTTMP/st_${L}_$$"
    if ! $BWA mem -t $THREADS -R "@RG\tID:${L}\tSM:${L}\tLB:${L}\tPL:ILLUMINA" "$REF" \
            "$FQDIR/${L}_1.fastq.gz" "$FQDIR/${L}_2.fastq.gz" 2>>"$LOG" \
          | $SAMTOOLS sort -@ $THREADS -m $SORT_MEM -T "$ST" -o "$BAMDIR/${L}.sorted.bam" - 2>>"$LOG"; then
      log "[$L] FAIL bwa|sort"; rm -rf "$FQDIR" "$BAMDIR" "${ST}"*; return 1
    fi
    rm -rf "$FQDIR" "${ST}"*
    log "[$L] sorted=$(du -h "$BAMDIR/${L}.sorted.bam"|cut -f1); MarkDuplicates"
    $SAMTOOLS index -@ $THREADS "$BAMDIR/${L}.sorted.bam"
    if ! $GATK --java-options "-Xmx32g" MarkDuplicates \
            -I "$BAMDIR/${L}.sorted.bam" -O "$DEDUP" -M "$BAMDIR/${L}.metrics.txt" \
            --TMP_DIR "$SORTTMP" 2>>"$LOG"; then
      log "[$L] FAIL MarkDuplicates"; rm -rf "$BAMDIR"; return 1
    fi
    rm -f "$BAMDIR/${L}.sorted.bam" "$BAMDIR/${L}.sorted.bam.bai"
    $SAMTOOLS index -@ $THREADS "$DEDUP"
    # Sanity: primary, mapped, non-duplicate reads should be a large fraction of ENA read_count
    local MAPPED; MAPPED=$($SAMTOOLS view -c -F 0xD04 "$DEDUP" 2>>"$LOG")
    log "[$L] primary-mapped-nondup reads=$MAPPED (ENA read_count=${READ_COUNT:-?})"
    if [ "${READ_COUNT:-0}" -gt 0 ] 2>/dev/null && [ "$MAPPED" -lt $(( READ_COUNT * 50 / 100 )) ]; then
      log "[$L] FAIL alignment: mapped $MAPPED < 50% of read_count $READ_COUNT — truncated input?"; rm -rf "$BAMDIR"; return 1
    fi
  fi

  local NBED; NBED=$(wc -l < "$BED")
  log "[$L] mpileup (-d $MPILEUP_D -q20 -Q20 --no-BAQ) over $NBED CDS C positions"
  if ! $SAMTOOLS mpileup -f "$REF" -l "$BED" -q 20 -Q 20 -d $MPILEUP_D --no-BAQ \
          "$DEDUP" 2>>"$LOG" > "${FINAL}.tmp"; then
    log "[$L] FAIL mpileup"; rm -f "${FINAL}.tmp"; return 1
  fi
  local NMP MINOK; NMP=$(wc -l < "${FINAL}.tmp"); MINOK=$(( NBED * 90 / 100 ))
  if [ "$NMP" -lt "$MINOK" ]; then
    log "[$L] FAIL mpileup incomplete: $NMP lines < $MINOK (90% of $NBED) — truncated or low coverage"; rm -f "${FINAL}.tmp"; return 1
  fi
  mv "${FINAL}.tmp" "$FINAL"
  echo "$NMP" > "${FINAL}.done"
  log "[$L] DONE mpileup=$(du -h "$FINAL"|cut -f1) lines=$NMP/$NBED"
  if [ -n "$GCS_OUT" ]; then
    if gsutil -q cp "$FINAL" "${GCS_OUT}/${L}_samtools.mpileup" && gsutil -q cp "${FINAL}.done" "${GCS_OUT}/${L}_samtools.mpileup.done"; then
      log "[$L] uploaded to ${GCS_OUT} (durable)"
    else
      log "[$L] WARN GCS upload failed — local copy kept at $FINAL"
    fi
  fi
  # Reclaim space only for clones we downloaded this run; keep the two pre-existing originals.
  if [ "$DOWNLOADED" -eq 1 ]; then rm -rf "$BAMDIR"; fi
  return 0
}

log "==== Path A v9 START — ${#ORDER[@]} clones queued: ${ORDER[*]} ===="
log "host=$(hostname) free=$(free_gb)G threads=$THREADS sort_mem=$SORT_MEM mpileup_d=$MPILEUP_D"
for L in "${ORDER[@]}"; do
  if [ -z "${URL[$L]:-}" ]; then log "[$L] UNKNOWN label, skipping"; continue; fi
  process_clone "$L" || log "[$L] clone failed — continuing to next"
done
log "==== Path A v9 COMPLETE — mpileups in $GOLD ===="
