#!/bin/bash
# Lei ASCERTAINMENT-MATCHED composition null: covered-but-UNEDITED TpC-C/GpA-G reference sites in the
# Detect-seq NT-rep1 assay universe (SRR11855272.bam), per-chrom parallel. awk keeps cov>=10 & NO mismatch
# (unedited) C/G positions, subsampled (NR%25==0). Output chrom/pos/ref -> python annotates TpC + RFD + OR.
B=/home/shaharh_quris_ai/miniconda3/envs/apobec/bin
REF=/data/ref/hg19.fa
BAM=/data/detectseq/SRR11855272.bam
OUT=/data/detectseq/covnull
mkdir -p $OUT
do_chrom(){
  ch=$1
  $B/samtools mpileup -r $ch -f $REF -q20 -Q20 -d 200 $BAM 2>/dev/null \
   | awk 'BEGIN{OFS="\t"} ($3=="C"||$3=="c"||$3=="G"||$3=="g") && $4>=10 && $5 !~ /[ACGTacgt]/ && (NR%25==0){print $1,$2,toupper($3)}' \
   > $OUT/${ch}.tsv
  echo "done $ch $(wc -l < $OUT/${ch}.tsv)"
}
export -f do_chrom; export B REF BAM OUT
printf "chr1\nchr2\nchr3\nchr7\nchr12\nchr17\n" | xargs -P 6 -I@ bash -c 'do_chrom @'
cat $OUT/chr*.tsv > $OUT/all.tsv
echo "COVNULL_CALL_DONE $(wc -l < $OUT/all.tsv)"
