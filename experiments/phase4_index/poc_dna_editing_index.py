"""POC: DNA editing index for BE4 (rAPOBEC1) off-target, Doman 2020 cohort, hg19.

MEASUREMENT POC (model deferred). Translates the Levanon RNA editing-index principles
to DNA off-target:
  - NO per-site coverage threshold: low-coverage sites are kept and ACCUMULATE via pooling.
  - Pooled (read-summed) rate estimator, treated minus matched control.
  - Many small edits across the genome accumulate into measurable genic burden.

Data in hand: 4 BE4-treated replicates + Parent (untreated) control.
  TREATED = BE4_clone1, BE4_clone5, BE4_clone6, BE4_clone7
  CONTROL = Parent_WGS

KEY DESIGN CHANGE vs parse_multiclone_pileup.py:
  The old script took the INTERSECTION of positions present in ALL samples, which silently
  drops a site that has 0 reads in even one clone -> violates the low-coverage principle.
  Here we POOL: a site is kept if (sum of treated coverage >= 1) AND (control coverage >= 1).
  A site seen in 3/4 clones still counts. Low-coverage sites accumulate across clones.

Germline whitelist: drop a site iff the Parent control VAF > 0.05 (inherited SNP, not low coverage).

Outputs (experiments/phase4_index/poc_results/):
  per_site_labels.parquet   - site-level (chrom,pos,gene,motif,strand, pooled treated/control edited+tot)
                              = the Phase-2 model's SITE-LEVEL LABELS
  per_gene_dei.parquet      - gene-level pooled DEI = the GENE-LEVEL LABELS
  per_clone_gene_rates.parquet - per-clone per-gene TpC rate (reproducibility)
  burden_capture.csv        - cumulative edit-burden captured by top-K genes (precision@K)
  qc_summary.txt            - headline numbers + pre-registered gates
"""
import os, sys, time
import numpy as np
import pandas as pd

DATA = '/tmp/poc_dna'
MP = f'{DATA}/mpileups'
BED = f'{DATA}/bed'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'poc_results')
os.makedirs(OUT, exist_ok=True)

TREATED = ['BE4_clone1', 'BE4_clone5', 'BE4_clone6', 'BE4_clone7']
CONTROL = ['Parent_WGS']
PARENT_VAF_CAP = 0.05
MIN_POS_PER_GENE = 20
N_TARGET = 15_666_970

def log(*a):
    print(f'[{time.strftime("%H:%M:%S")}]', *a, file=sys.stderr, flush=True)

# ---------------------------------------------------------------------------
# 1. Load BED annotation -> position universe (chrom,pos) -> index, gene/motif/strand
# ---------------------------------------------------------------------------
log('loading BED annotation...')
pos_index = {}          # (chrom,pos1) -> idx
genes = []; motifs = []; strands = []
for bedp, mot in [(f'{BED}/cds_tpc.bed', 'tpc'), (f'{BED}/cds_npc.bed', 'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            chrom, s, e, gene, sc, strand = f
            pos1 = int(s) + 1
            key = (chrom, pos1)
            if key in pos_index:
                continue
            pos_index[key] = len(genes)
            genes.append(gene); motifs.append(mot); strands.append(strand)
N = len(genes)
genes = np.array(genes, dtype=object)
motifs = np.array(motifs)
strands = np.array(strands)
is_tpc = (motifs == 'tpc')
strand_plus = (strands == '+')
log(f'annotated positions: {N:,}  (tpc={is_tpc.sum():,}  npc={(~is_tpc).sum():,})')


def parse_bases_counts(bases, ref):
    """Return (A,C,G,T) counts from an mpileup read-bases string."""
    a = c = g = t = 0
    R = ref.upper()
    i = 0; n = len(bases)
    while i < n:
        ch = bases[i]
        if ch == '^':
            i += 2; continue
        if ch == '$':
            i += 1; continue
        if ch in '+-':
            j = i + 1
            while j < n and bases[j].isdigit():
                j += 1
            try:
                i = j + int(bases[i+1:j])
            except ValueError:
                i += 1
            continue
        if ch == '*':
            i += 1; continue
        if ch in '.,':
            if R == 'A': a += 1
            elif R == 'C': c += 1
            elif R == 'G': g += 1
            elif R == 'T': t += 1
        else:
            u = ch.upper()
            if u == 'A': a += 1
            elif u == 'C': c += 1
            elif u == 'G': g += 1
            elif u == 'T': t += 1
        i += 1
    return a, c, g, t


def parse_sample(label):
    """Fill arrays (len N) of total, edited (C->T), noise (other mismatch) for this sample.

    + strand ref=C: edited=T,            noise=(A+G)/2
    - strand ref=G: edited=A (==C->T),   noise=(C+T)/2
    Positions not in the BED universe, or rows whose ref doesn't match the strand, are skipped.
    """
    tot = np.zeros(N, dtype=np.int32)
    edited = np.zeros(N, dtype=np.int32)
    noise = np.zeros(N, dtype=np.float32)
    path = f'{MP}/{label}_samtools.mpileup'
    nrows = 0; kept = 0
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            nrows += 1
            key = (f[0], int(f[1]))
            idx = pos_index.get(key)
            if idx is None:
                continue
            ref = f[2].upper()
            a, c, g, t = parse_bases_counts(f[4], ref)
            tt = a + c + g + t
            if tt < 1:
                continue
            if strand_plus[idx]:
                if ref != 'C':
                    continue
                ed = t; ns = (a + g) / 2.0
            else:
                if ref != 'G':
                    continue
                ed = a; ns = (c + t) / 2.0
            tot[idx] = tt; edited[idx] = ed; noise[idx] = ns
            kept += 1
            if nrows % 3_000_000 == 0:
                log(f'  {label}: {nrows:,} rows -> kept {kept:,}')
    log(f'  {label}: total {nrows:,} kept {kept:,}')
    return tot, edited, noise


# ---------------------------------------------------------------------------
# 2. Parse all samples into aligned arrays
# ---------------------------------------------------------------------------
S = {}
for s in TREATED + CONTROL:
    log(f'parsing {s}...')
    S[s] = parse_sample(s)

# Pooled treated, control
t_tot = np.zeros(N, dtype=np.int64)
t_ed = np.zeros(N, dtype=np.int64)
t_noise = np.zeros(N, dtype=np.float64)
for s in TREATED:
    tot, ed, ns = S[s]
    t_tot += tot; t_ed += ed; t_noise += ns
c_tot = S['Parent_WGS'][0].astype(np.int64)
c_ed = S['Parent_WGS'][1].astype(np.int64)
c_noise = S['Parent_WGS'][2].astype(np.float64)

# ---------------------------------------------------------------------------
# 3. POOLED keep mask (Levanon-faithful: no per-site coverage floor beyond >=1 pooled)
#    + germline whitelist (Parent VAF > cap dropped)
# ---------------------------------------------------------------------------
covered = (t_tot >= 1) & (c_tot >= 1)
parent_vaf = np.divide(c_ed, c_tot, out=np.zeros(N), where=c_tot > 0)
germline = parent_vaf > PARENT_VAF_CAP
keep = covered & (~germline)
log(f'covered (pooled t>=1 & c>=1): {covered.sum():,} / {N:,}')
log(f'germline dropped (Parent VAF>{PARENT_VAF_CAP}): {germline[covered].sum():,}')
log(f'kept sites: {keep.sum():,}')

# ---------------------------------------------------------------------------
# 4. Per-site label table (the model's site-level labels)
# ---------------------------------------------------------------------------
chroms = np.empty(N, dtype=object)
poss = np.empty(N, dtype=np.int64)
for (chrom, p1), idx in pos_index.items():
    chroms[idx] = chrom; poss[idx] = p1
site = pd.DataFrame({
    'chrom': chroms[keep], 'pos': poss[keep], 'gene': genes[keep],
    'motif': motifs[keep], 'strand': strands[keep],
    't_edited': t_ed[keep], 't_tot': t_tot[keep],
    'c_edited': c_ed[keep], 'c_tot': c_tot[keep],
})
site['t_rate'] = site.t_edited / site.t_tot
site['c_rate'] = site.c_edited / site.c_tot
site['site_dei'] = site.t_rate - site.c_rate
site.to_parquet(f'{OUT}/per_site_labels.parquet')
log(f'wrote per_site_labels.parquet: {len(site):,} sites')

# ---------------------------------------------------------------------------
# 5. Per-gene pooled DEI
# ---------------------------------------------------------------------------
df = pd.DataFrame({
    'gene': genes[keep], 'is_tpc': is_tpc[keep],
    't_ed': t_ed[keep], 't_tot': t_tot[keep], 't_noise': t_noise[keep],
    'c_ed': c_ed[keep], 'c_tot': c_tot[keep], 'c_noise': c_noise[keep],
})
def agg_gene(sub):
    out = {}
    for mo, tag in [(True, 'tpc'), (False, 'npc')]:
        m = sub.is_tpc.values == mo
        tt = sub.t_tot.values[m].sum(); te = sub.t_ed.values[m].sum(); tn = sub.t_noise.values[m].sum()
        ct = sub.c_tot.values[m].sum(); ce = sub.c_ed.values[m].sum(); cn = sub.c_noise.values[m].sum()
        out[f'{tag}_npos'] = int(m.sum())
        out[f'{tag}_t_rate'] = te / tt if tt else 0.0
        out[f'{tag}_c_rate'] = ce / ct if ct else 0.0
        out[f'{tag}_t_noise'] = tn / tt if tt else 0.0
        out[f'{tag}_t_excess_edits'] = te - (ce / ct * tt if ct else 0.0)  # excess edits over control rate
    return pd.Series(out)
gene_df = df.groupby('gene').apply(agg_gene).reset_index()
gene_df['DEI_tpc'] = gene_df.tpc_t_rate - gene_df.tpc_c_rate
gene_df['DEI_npc'] = gene_df.npc_t_rate - gene_df.npc_c_rate
gene_df['noise_gap_tpc'] = gene_df.tpc_t_rate - gene_df.tpc_t_noise
gene_full = gene_df.copy()
gene_df = gene_df[(gene_df.tpc_npos >= MIN_POS_PER_GENE)].copy()
gene_df.to_parquet(f'{OUT}/per_gene_dei.parquet')
log(f'wrote per_gene_dei.parquet: {len(gene_df):,} genes (>= {MIN_POS_PER_GENE} TpC pos)')

# ---------------------------------------------------------------------------
# 6. Reproducibility: per-clone per-gene TpC rate
# ---------------------------------------------------------------------------
rep = {}
gidx_tpc = is_tpc & keep
for s in TREATED:
    tot, ed, _ = S[s]
    sub = pd.DataFrame({'gene': genes[gidx_tpc], 'ed': ed[gidx_tpc], 'tot': tot[gidx_tpc]})
    gr = sub.groupby('gene').sum()
    rep[s] = (gr.ed / gr.tot.replace(0, np.nan))
rep_df = pd.DataFrame(rep)
rep_df = rep_df.loc[gene_df.gene.values].dropna()
rep_df.to_parquet(f'{OUT}/per_clone_gene_rates.parquet')
corr = rep_df.corr(method='spearman')
log('per-clone per-gene TpC-rate Spearman correlation:')
log('\n' + corr.round(3).to_string())

# ---------------------------------------------------------------------------
# 7. Burden-capture curve (precision@K): cumulative excess-edit burden by top-K genes
# ---------------------------------------------------------------------------
bc = gene_df[['gene', 'DEI_tpc', 'tpc_t_excess_edits', 'tpc_npos']].copy()
bc = bc.sort_values('tpc_t_excess_edits', ascending=False).reset_index(drop=True)
total_excess = bc.tpc_t_excess_edits.clip(lower=0).sum()
bc['cum_excess'] = bc.tpc_t_excess_edits.clip(lower=0).cumsum()
bc['cum_frac_burden'] = bc.cum_excess / total_excess
bc['rank'] = np.arange(1, len(bc) + 1)
bc['frac_genes'] = bc['rank'] / len(bc)
bc[['rank', 'gene', 'DEI_tpc', 'tpc_t_excess_edits', 'cum_frac_burden', 'frac_genes']].to_csv(
    f'{OUT}/burden_capture.csv', index=False)


def frac_genes_for(target):
    hit = bc[bc.cum_frac_burden >= target]
    return (hit['rank'].iloc[0], hit['frac_genes'].iloc[0]) if len(hit) else (len(bc), 1.0)

# ---------------------------------------------------------------------------
# 8. QC summary + pre-registered gates
# ---------------------------------------------------------------------------
def stat(a):
    a = np.asarray(a, dtype=float)
    return f'n={len(a)} mean={a.mean():.4g} med={np.median(a):.4g} q90={np.quantile(a,0.9):.4g}'

g = gene_df
lines = [
    '=== POC: BE4 DNA editing index (4 BE4 clones vs Parent, hg19) ===',
    f'treated={TREATED}',
    f'control={CONTROL}',
    f'estimator=POOLED (Levanon: no per-site coverage floor; low-cov sites accumulate)',
    f'PARENT_VAF_CAP={PARENT_VAF_CAP}  MIN_POS_PER_GENE={MIN_POS_PER_GENE}',
    '',
    f'sites: covered={covered.sum():,}  kept(after germline)={keep.sum():,}  genes(>= {MIN_POS_PER_GENE} TpC pos)={len(g):,}',
    '',
    '--- pooled per-gene rates ---',
    f'Treated TpC rate:    {stat(g.tpc_t_rate)}',
    f'Control TpC rate:    {stat(g.tpc_c_rate)}',
    f'Treated nonTpC rate: {stat(g.npc_t_rate)}',
    f'Control nonTpC rate: {stat(g.npc_c_rate)}',
    '',
    '--- DEI (HEADLINE) ---',
    f'DEI_tpc: {stat(g.DEI_tpc)}',
    f'DEI_npc (motif-negative control): {stat(g.DEI_npc)}',
    '',
    '--- reproducibility (per-clone per-gene TpC rate, Spearman) ---',
    corr.round(3).to_string(),
    f'mean off-diagonal Spearman: {corr.values[np.triu_indices_from(corr.values,1)].mean():.3f}',
    '',
    '--- pre-registered gates ---',
    f'[a] DEI_tpc mean > 0:                 {g.DEI_tpc.mean() > 0}  ({g.DEI_tpc.mean():.3g})',
    f'[b] mean(DEI_tpc)/mean(DEI_npc):       {g.DEI_tpc.mean()/max(g.DEI_npc.mean(),1e-12):.2f}x  (>2 convincing, >1.3 weak)',
    f'[c] frac genes DEI_tpc > DEI_npc:      {(g.DEI_tpc > g.DEI_npc).mean():.3f}  (>0.70 target)',
    f'[d] frac genes treated TpC > control:  {(g.tpc_t_rate > g.tpc_c_rate).mean():.3f}',
    '',
    '--- top-100 motif specificity ---',
    f'top-100 DEI_tpc mean: {g.nlargest(100,"DEI_tpc").DEI_tpc.mean():.4g}',
    f'top-100 DEI_npc mean: {g.nlargest(100,"DEI_tpc").DEI_npc.mean():.4g}',
    f'top-100 ratio: {g.nlargest(100,"DEI_tpc").DEI_tpc.mean()/max(g.nlargest(100,"DEI_tpc").DEI_npc.mean(),1e-12):.2f}x',
    '',
    '--- burden-capture (precision@K) ---',
]
for tgt in (0.50, 0.80, 0.90):
    rnk, fg = frac_genes_for(tgt)
    lines.append(f'  {int(tgt*100)}% of excess TpC edit burden captured by top {rnk:,} genes ({100*fg:.1f}% of genes)')

text = '\n'.join(str(x) for x in lines)
print(text)
with open(f'{OUT}/qc_summary.txt', 'w') as fh:
    fh.write(text + '\n')
log('wrote qc_summary.txt')
