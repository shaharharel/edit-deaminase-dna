#!/usr/bin/env bash
# Run on ai-chem. Computes the DNA editing index for Doman BE4_clone1 vs Parent_WGS.
# Expects samtools in PATH and refs at /mnt/data/ref/hg38/hg38.fa + /mnt/data/ref/refGene.txt.
set -euo pipefail

REF=/mnt/data/ref/hg38/hg38.fa
REFGENE=/mnt/data/ref/refGene.txt
OUTDIR=/mnt/data/dna_features/index
mkdir -p "$OUTDIR"
NPROC=${NPROC:-32}

# 1) Build TpC + non-TpC CDS BED (motif-negative control set is the non-TpC C's in CDS)
python3 scripts/data/build_tpc_cds_bed.py "$REF" "$REFGENE" "$OUTDIR"

# 2) mpileup at TpC + non-TpC positions, for BE4_clone1 and Parent_WGS
#    -q 20 -Q 20  : map and base quality filters
#    -d 200       : max depth per position (avoid pile-up explosion at duplications)
#    --no-BAQ     : don't do BAQ realignment (faster; we don't need indel realign)
cat "$OUTDIR/cds_tpc.bed" "$OUTDIR/cds_npc.bed" | cut -f1-3 | sort -k1,1 -k2,2n -u > "$OUTDIR/positions.bed"
echo "positions to pileup: $(wc -l < "$OUTDIR/positions.bed")"

for SAMPLE in BE4_clone1 Parent_WGS; do
  BAM=/mnt/data/doman_2020/${SAMPLE}.bam
  if [[ ! -f $BAM ]]; then echo "MISSING $BAM" >&2; continue; fi
  echo "===> mpileup $SAMPLE"
  samtools mpileup -f "$REF" -l "$OUTDIR/positions.bed" -q 20 -Q 20 -d 200 --no-BAQ "$BAM" \
    > "$OUTDIR/${SAMPLE}.mpileup" 2> "$OUTDIR/${SAMPLE}.mpileup.log"
  echo "  wrote $OUTDIR/${SAMPLE}.mpileup ($(wc -l < "$OUTDIR/${SAMPLE}.mpileup") rows)"
done

# 3) Parse mpileup into per-position counts (with mismatch-direction control)
for SAMPLE in BE4_clone1 Parent_WGS; do
  echo "===> parse $SAMPLE"
  python3 scripts/data/parse_pileup_for_index.py \
    "$OUTDIR/${SAMPLE}.mpileup" "$OUTDIR/cds_tpc.bed" "$OUTDIR/cds_npc.bed" \
    "$OUTDIR/${SAMPLE}.counts.parquet" "$SAMPLE"
done

# 4) Compute per-gene editing index (low-VAF filter + mismatch-direction + treated-control)
echo "===> compute DEI"
python3 experiments/phase4_index/compute_dna_editing_index.py \
  "$OUTDIR/BE4_clone1.counts.parquet" "$OUTDIR/Parent_WGS.counts.parquet" \
  "$OUTDIR/per_gene_index.parquet" "$OUTDIR/qc_summary.txt" \
  0.40 20

# 5) Compute per-gene f x g prediction
echo "===> per-gene f x g prediction"
python3 experiments/phase4_index/predict_per_gene_fxg.py \
  "$OUTDIR/cds_tpc.bed" "$REF" "/mnt/data/dna_features/bins_1mb_v3.parquet" \
  "/mnt/data/dna_features/canonical_sites.parquet" \
  "$OUTDIR/per_gene_fxg.parquet"

# 6) Validate empirical DEI vs predicted f x g
echo "===> validate"
python3 experiments/phase4_index/validate_index_vs_fxg.py \
  "$OUTDIR/per_gene_index.parquet" "$OUTDIR/per_gene_fxg.parquet" \
  "$OUTDIR/per_gene_index_vs_fxg.parquet" | tee "$OUTDIR/validation.txt"

echo
echo "ALL DONE. Outputs in $OUTDIR/"
echo "  per_gene_index.parquet         - Levanon-style DEI per gene (empirical)"
echo "  per_gene_fxg.parquet           - f x g predicted per-gene score"
echo "  per_gene_index_vs_fxg.parquet  - merged + spearman/recall"
echo "  qc_summary.txt                  - motif-spectrum + VAF-filter sanity"
echo "  validation.txt                  - predicted vs empirical correlation"
