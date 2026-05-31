"""Parse samtools mpileup -> per-position read counts at editable C/G positions, with
mismatch-direction control. Inputs: pileup file from BAM, bed of TpC and bed of nonTpC.
Output: per-position parquet with cols
  chrom, pos, strand, motif (tpc|npc), gene, depth,
  n_C, n_T, n_A, n_G  (read counts on the forward strand; strand-corrected)
  vaf_main         (variant allele freq of the EDITED direction: C->T on '+' strand, G->A on '-')
  vaf_noise        (mismatch-direction control: G->A on '+' strand, C->T on '-' strand)

Use it for per-gene aggregation downstream. Designed to handle pileup of ~9.7M sites x 2 BAMs.
"""
import sys, pandas as pd, numpy as np, re
from collections import defaultdict

PILEUP = sys.argv[1]            # samtools mpileup output (one row per requested position)
TPC_BED = sys.argv[2]           # cds_tpc.bed (chrom,start,end,gene,score,strand)
NPC_BED = sys.argv[3]           # cds_npc.bed
OUT_PARQUET = sys.argv[4]       # output parquet path
SAMPLE = sys.argv[5] if len(sys.argv) > 5 else 'sample'

# load motif & gene annotation for each position (chrom, pos1based, strand) -> (gene, motif_class)
ann = {}
for bedp, mot in [(TPC_BED,'tpc'),(NPC_BED,'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom, s, e, gene, sc, strand = f
            ann[(chrom, int(s)+1, strand)] = (gene, mot)

# pileup parsing: count A C G T from the bases string (5th column)
# samtools format: chrom pos ref depth bases quals
# bases use '.' = match-fwd-strand, ',' = match-rev-strand, ACGT = mismatch-fwd, acgt = mismatch-rev
# we want strand-aware mismatch counts. We pull motif strand from ann; we count reads independently of read direction
# (a mismatch to ref on either read direction counts as the variant).

def parse_bases(bases, ref):
    # collapse insertions/deletions to ignore them, drop ^X start-of-read marker (and its mapq),
    # drop '$' end-of-read marker
    out = {'A':0,'C':0,'G':0,'T':0}
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == '^':  # start of read followed by mapq char
            i += 2; continue
        if c == '$':
            i += 1; continue
        if c in '+-':  # indel length follows: e.g. +2AC
            j = i+1
            while j < len(bases) and bases[j].isdigit():
                j += 1
            ln = int(bases[i+1:j])
            i = j + ln
            continue
        if c in '.,':
            out[ref.upper()] = out.get(ref.upper(),0) + 1
        elif c.upper() in 'ACGT':
            out[c.upper()] = out.get(c.upper(),0) + 1
        i += 1
    return out

rows = []
with open(PILEUP) as fh:
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 5: continue
        chrom, pos, ref, depth, bases = f[0], int(f[1]), f[2], int(f[3]), f[4]
        # try both strands (a position may be in ann under '+' or '-')
        for strand in ('+','-'):
            key = (chrom, pos, strand)
            if key not in ann: continue
            gene, motif = ann[key]
            counts = parse_bases(bases, ref)
            n_A = counts.get('A',0); n_C = counts.get('C',0); n_G = counts.get('G',0); n_T = counts.get('T',0)
            tot = n_A + n_C + n_G + n_T
            if tot == 0: continue
            # for '+' strand TpC: edited direction = C->T, noise = G->A
            # for '-' strand TpC (i.e. GpA on fwd): edited direction (on fwd) = G->A, noise = C->T
            if strand == '+':
                v_main = n_T / tot
                v_noise = n_A / tot   # A at a C position; A->C would be the symmetric noise here
                # Actually for C->T direction, the natural noise floor is the OTHER 2 mismatches:
                # use the avg of the two non-edited mismatches (A and G at a C ref)
                v_noise = (n_A + n_G) / (2 * tot)
            else:
                v_main = n_A / tot    # G->A on fwd strand = TpC -> TpT on the reverse
                v_noise = (n_C + n_T) / (2 * tot)
            rows.append({
                'chrom':chrom, 'pos':pos, 'strand':strand, 'motif':motif, 'gene':gene,
                'depth':tot, 'n_A':n_A, 'n_C':n_C, 'n_G':n_G, 'n_T':n_T,
                'vaf_main':v_main, 'vaf_noise':v_noise, 'sample':SAMPLE
            })

df = pd.DataFrame(rows)
df.to_parquet(OUT_PARQUET, compression='snappy')
print(f"{SAMPLE}: {len(df):,} rows | TpC {(df.motif=='tpc').sum():,} | nonTpC {(df.motif=='npc').sum():,}", file=sys.stderr)
print(f"  median depth: {df.depth.median():.1f}", file=sys.stderr)
print(f"  mean VAF_main (all): {df.vaf_main.mean():.4f}", file=sys.stderr)
print(f"  mean VAF_main (low-VAF<0.4): {df[df.vaf_main<0.4].vaf_main.mean():.5f}", file=sys.stderr)
