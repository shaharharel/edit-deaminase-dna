"""Build BED of TpC dinucleotide positions in CDS (both strands).
Each row = a single C that is the edited base of a TpC motif on its strand.
Output: cds_tpc.bed (chrom, pos0, pos1, gene_name, score, strand)
Used as the mpileup -l target list for the DNA editing index.

QC: also writes cds_npc.bed (non-TpC C's in CDS) for the motif-negative control,
and per-gene TpC count summary cds_tpc_per_gene.tsv.
"""
import pandas as pd, numpy as np
from pyfaidx import Fasta
from collections import defaultdict
import sys, os

HG38 = sys.argv[1] if len(sys.argv) > 1 else '/mnt/data/ref/hg38/hg38.fa'
REFGENE = sys.argv[2] if len(sys.argv) > 2 else '/mnt/data/ref/refGene.txt'
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else '/mnt/data/dna_features'

fa = Fasta(HG38)
CH = [f'chr{i}' for i in range(1,23)] + ['chrX']

# refGene: parse CDS exon coords per gene
# columns: 0=bin 1=name 2=chrom 3=strand 4=txStart 5=txEnd 6=cdsStart 7=cdsEnd
# 8=exonCount 9=exonStarts 10=exonEnds 12=name2 (gene symbol)
rg = pd.read_csv(REFGENE, sep='\t', header=None,
    names=['bn','nm','chrom','strand','ts','te','cs','ce','ec','es','ee','sc','name2','a','c','f'])
rg = rg[(rg.chrom.isin(CH)) & (rg.cs < rg.ce)].copy()
print(f"refGene rows in CH: {len(rg)}", file=sys.stderr)

# Build the union CDS-intersect-exon intervals per gene
def cds_exons(r):
    starts = list(map(int, r.es.rstrip(',').split(',')))
    ends   = list(map(int, r.ee.rstrip(',').split(',')))
    out = []
    for s,e in zip(starts, ends):
        a = max(s, r.cs); b = min(e, r.ce)
        if a < b: out.append((a,b))
    return out

tpc_rows = []
npc_rows = []
per_gene_tpc = defaultdict(int)
per_gene_npc = defaultdict(int)

# PERF FIX: load each chromosome ONCE, sort refGene rows by chrom so we
# don't re-load the chromosome FASTA per transcript (was loading hundreds
# of times — O(N_tx * chrom_len) ~ 70k * 200M for chr1 alone, hours).
rg = rg.sort_values(['chrom','ts']).reset_index(drop=True)
cur_chrom = None
seq = None
for r in rg.itertuples():
    if r.chrom != cur_chrom:
        seq = fa[r.chrom][:].seq.upper()
        cur_chrom = r.chrom
        print(f"  loaded {r.chrom} ({len(seq):,} bp)", file=sys.stderr)
    for a,b in cds_exons(r):
        # QA-H9 fix: index into full chromosome seq so we can check upstream/downstream
        # bases at exon edges (and exclude CpG from npc since CpG spontaneous deamination
        # contaminates the motif-negative control).
        for pos in range(max(a,1), min(b, len(seq)-1)):
            base = seq[pos]
            up = seq[pos-1]
            dn = seq[pos+1]
            if base == 'C':
                # CpG on + strand: C with G downstream -> exclude (spontaneous deamination)
                if dn == 'G':
                    continue
                if up == 'T':
                    tpc_rows.append((r.chrom, pos, pos+1, r.name2, 0, '+'))
                    per_gene_tpc[(r.chrom, r.name2)] += 1
                else:
                    npc_rows.append((r.chrom, pos, pos+1, r.name2, 0, '+'))
                    per_gene_npc[(r.chrom, r.name2)] += 1
            elif base == 'G':
                # The C on the - strand at this position; CpG on - strand means upstream is C on + strand
                # (i.e. + strand has CG palindrome at pos-1,pos). Exclude.
                if up == 'C':
                    continue
                # TpC on the reverse strand = GpA on this strand
                if dn == 'A':
                    tpc_rows.append((r.chrom, pos, pos+1, r.name2, 0, '-'))
                    per_gene_tpc[(r.chrom, r.name2)] += 1
                else:
                    npc_rows.append((r.chrom, pos, pos+1, r.name2, 0, '-'))
                    per_gene_npc[(r.chrom, r.name2)] += 1

os.makedirs(OUTDIR, exist_ok=True)
tpc_df = pd.DataFrame(tpc_rows, columns=['chrom','start','end','gene','score','strand']).drop_duplicates(['chrom','start','strand'])
npc_df = pd.DataFrame(npc_rows, columns=['chrom','start','end','gene','score','strand']).drop_duplicates(['chrom','start','strand'])

tpc_df.to_csv(f'{OUTDIR}/cds_tpc.bed', sep='\t', header=False, index=False)
npc_df.to_csv(f'{OUTDIR}/cds_npc.bed', sep='\t', header=False, index=False)

per_gene = pd.DataFrame([
    {'chrom':c, 'gene':g, 'n_tpc': per_gene_tpc[(c,g)], 'n_npc': per_gene_npc[(c,g)]}
    for (c,g) in set(list(per_gene_tpc.keys()) + list(per_gene_npc.keys()))
])
per_gene.to_csv(f'{OUTDIR}/cds_tpc_per_gene.tsv', sep='\t', index=False)

print(f"TpC sites: {len(tpc_df):>10,}  ({(tpc_df.strand=='+').sum():,} + / {(tpc_df.strand=='-').sum():,} -)", file=sys.stderr)
print(f"non-TpC C's: {len(npc_df):>10,}  (motif-negative control set)", file=sys.stderr)
print(f"genes: {per_gene.gene.nunique():,}", file=sys.stderr)
print(f"median TpC per gene: {per_gene.n_tpc.median():.0f}, p10={per_gene.n_tpc.quantile(0.10):.0f}, p90={per_gene.n_tpc.quantile(0.90):.0f}", file=sys.stderr)
