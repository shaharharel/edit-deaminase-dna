"""24-clone cross-deaminase DNA editing-index POC (Doman 2020, hg19) — FULL cohort.

Successor to poc_12clone.py: now the complete Doman cross-deaminase WGS set.
Cohorts (24 clones, all durable in gs://ai-temp/apobec-genome-cache/mpileups):
  BE4      = clone1..8       (rAPOBEC1, dirty)         n=8
  YE1-BE4  = clone1..8       (engineered safer)        n=8
  nCas9    = clone1..7       (deaminase-free control)  n=7
  Parent   = Parent_WGS      (untreated baseline)      n=1

Levanon-faithful: pooled (read-summed) RATIO estimator, NO per-site coverage floor,
germline whitelist via Parent VAF>0.05. Motif channels TpC (APOBEC1 sig) vs non-TpC
(per-gene internal negative control). DEI = treated_ratio - parent_ratio.

Same analyses A-G as poc_12clone.py; split-half within BE4 is now 4v4 ({1,2,3,4} vs {5,6,7,8}).
Reads /tmp/poc_dna/{mpileups,bed}. Writes experiments/phase4_index/poc24_results/.
Per-site label table retains chrom,pos,strand,gene,motif (+ per-cohort ed/tot) so Evo
sequence-context windows and structure features can bolt on as extra columns for the
Phase-2 model (Architecture B).
"""
import os, sys, time
import numpy as np, pandas as pd

DATA = '/tmp/poc_dna'; MP = f'{DATA}/mpileups'; BED = f'{DATA}/bed'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'poc24_results')
os.makedirs(OUT, exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)

COHORTS = {
    'BE4':     [f'BE4_clone{i}' for i in range(1, 9)],
    'YE1-BE4': [f'YE1-BE4_clone{i}' for i in range(1, 9)],
    'nCas9':   [f'nCas9_clone{i}' for i in range(1, 8)],
}
PARENT = 'Parent_WGS'
ALL_CLONES = [c for v in COHORTS.values() for c in v] + [PARENT]
PARENT_VAF_CAP = 0.05
MIN_POS = 20

def log(*a): print(f'[{time.strftime("%H:%M:%S")}]', *a, file=sys.stderr, flush=True)

# ---------------- load BED -> position index ----------------
log('loading BED...')
pos_index = {}; chroms = []; positions = []; genes = []; motifs = []; strands = []
for bedp, mot in [(f'{BED}/cds_tpc.bed', 'tpc'), (f'{BED}/cds_npc.bed', 'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom = f[0]; pos = int(f[1]) + 1
            key = (chrom, pos)
            if key in pos_index: continue
            pos_index[key] = len(genes)
            chroms.append(chrom); positions.append(pos); genes.append(f[3]); motifs.append(mot); strands.append(f[5])
N = len(genes)
chroms = np.array(chroms, dtype=object); positions = np.array(positions, dtype=np.int32)
genes = np.array(genes, dtype=object); is_tpc = (np.array(motifs) == 'tpc'); strand_plus = (np.array(strands) == '+')
log(f'positions N={N:,} (tpc={is_tpc.sum():,} npc={(~is_tpc).sum():,})')

# ---------------- mpileup parse ----------------
def parse_counts(bases, ref):
    a = c = g = t = 0; R = ref.upper(); i = 0; n = len(bases)
    while i < n:
        ch = bases[i]
        if ch == '^': i += 2; continue
        if ch == '$': i += 1; continue
        if ch in '+-':
            j = i + 1
            while j < n and bases[j].isdigit(): j += 1
            try: i = j + int(bases[i+1:j])
            except ValueError: i += 1
            continue
        if ch == '*': i += 1; continue
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
    """(tot, ed) int32 over N. edit = C->T equiv: +strand refC->T, -strand refG->A."""
    tot = np.zeros(N, dtype=np.int32); ed = np.zeros(N, dtype=np.int32)
    path = f'{MP}/{label}_samtools.mpileup'
    nlines = 0
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 5: continue
            idx = pos_index.get((f[0], int(f[1])))
            if idx is None: continue
            ref = f[2].upper(); a, c, g, t = parse_counts(f[4], ref); tt = a + c + g + t
            if tt < 1: continue
            if strand_plus[idx]:
                if ref != 'C': continue
                e = t
            else:
                if ref != 'G': continue
                e = a
            tot[idx] = tt; ed[idx] = e; nlines += 1
    log(f'  parsed {label}: {nlines:,} usable positions')
    return tot, ed

log('parsing all 24 clones (slow part)...')
CACHE = f'{DATA}/parse_cache_24clone.npz'
S = {}
if os.path.exists(CACHE):
    log(f'loading parse cache {CACHE} ...')
    z = np.load(CACHE, allow_pickle=False)
    if int(z['N']) == N and all(f'{cl}__tot' in z.files for cl in ALL_CLONES):
        for cl in ALL_CLONES: S[cl] = (z[f'{cl}__tot'], z[f'{cl}__ed'])
        log('parse cache OK -> skipped parse')
    else:
        log('parse cache stale; reparsing'); S = {}
if not S:
    # parse all clones in parallel (fork workers inherit pos_index/strand_plus/N) —
    # uses the box's cores to shrink the single-threaded window (spot-preemption safety)
    from multiprocessing import Pool
    # CAP at 6: each fork worker touches the ~2.5GB shared pos_index dict (refcount
    # COW), so concurrency * 2.5GB must stay well under RAM. 24-way wedged a 128GB box.
    nproc = min(6, len(ALL_CLONES), max(1, (os.cpu_count() or 4)))
    log(f'parsing {len(ALL_CLONES)} clones across {nproc} procs (COW-capped)...')
    with Pool(nproc, maxtasksperchild=1) as pool:
        results = pool.map(parse_sample, ALL_CLONES, chunksize=1)
    S = dict(zip(ALL_CLONES, results))
    np.savez(CACHE, N=np.int64(N),
             **{f'{cl}__tot': S[cl][0] for cl in ALL_CLONES},
             **{f'{cl}__ed': S[cl][1] for cl in ALL_CLONES})
    log(f'saved parse cache -> {CACHE}')

ptot, ped = S[PARENT]
parent_vaf = np.divide(ped, ptot, out=np.zeros(N), where=ptot > 0)
keep = (ptot >= 1) & (parent_vaf <= PARENT_VAF_CAP)
log(f'kept after germline filter: {keep.sum():,}/{N:,}')

def pool(clones):
    tot = np.zeros(N, dtype=np.int64); ed = np.zeros(N, dtype=np.int64)
    for cl in clones: tot += S[cl][0]; ed += S[cl][1]
    return tot, ed

def bulk_rate(clones, mask):
    tot, ed = pool(clones); m = keep & mask
    T = tot[m].sum(); E = ed[m].sum()
    return (E / T if T else 0.0), int(E), int(T)

report = []
report.append("=== 24-CLONE CROSS-DEAMINASE DNA EDITING-INDEX POC (full Doman cohort) ===")
report.append(f"positions N={N:,}; kept(germline)={int(keep.sum()):,}")
report.append("cohorts: BE4 n=8, YE1-BE4 n=8, nCas9 n=7, Parent n=1")
report.append("\n--- B. BULK pooled C->T rate (TpC channel), per cohort ---")
prate, pe, pt = bulk_rate([PARENT], is_tpc)
report.append(f"  Parent : {prate:.3e}  ({pe:,}/{pt:,})")
for name, clones in COHORTS.items():
    r, e, t = bulk_rate(clones, is_tpc)
    report.append(f"  {name:8s}: {r:.3e}  ({e:,}/{t:,})  vs Parent x{r/prate:.3f}")

report.append("\n--- D. MOTIF SPECIFICITY: pooled TpC rate / nonTpC rate (>1 = APOBEC-like) ---")
for name, clones in ([('Parent', [PARENT])] + list(COHORTS.items())):
    rt, _, _ = bulk_rate(clones, is_tpc); rn, _, _ = bulk_rate(clones, ~is_tpc)
    spec = rt / rn if rn else float('nan')
    report.append(f"  {name:8s}: TpC={rt:.3e} nonTpC={rn:.3e}  specificity={spec:.3f}  (TpC-nonTpC={rt-rn:+.3e})")

def per_gene_dei(clones):
    tot, ed = pool(clones); k = keep & (tot >= 1)
    df = pd.DataFrame({'gene': genes[k], 'tpc': is_tpc[k], 'ed': ed[k], 'tot': tot[k],
                       'c_ed': ped[k], 'c_tot': ptot[k]})
    rows = []
    for gene, sub in df.groupby('gene'):
        m = sub.tpc.values
        for tag, sel in [('tpc', m), ('npc', ~m)]:
            tt = sub.tot.values[sel].sum(); te = sub.ed.values[sel].sum()
            ct = sub.c_tot.values[sel].sum(); ce = sub.c_ed.values[sel].sum()
            rows.append((gene, tag, int(sel.sum()), te / tt if tt else 0.0, ce / ct if ct else 0.0))
    r = pd.DataFrame(rows, columns=['gene', 'tag', 'npos', 't_rate', 'c_rate'])
    r['dei'] = r.t_rate - r.c_rate
    piv = r.pivot(index='gene', columns='tag', values='dei')
    piv['npos_tpc'] = r[r.tag == 'tpc'].set_index('gene').npos
    piv['mc'] = piv['tpc'] - piv['npc']
    return piv

log('per-gene DEI per cohort...')
gene_tables = {}
for name, clones in COHORTS.items():
    piv = per_gene_dei(clones); piv = piv[piv.npos_tpc >= MIN_POS]
    gene_tables[name] = piv; piv.to_parquet(f'{OUT}/per_gene_dei_{name}.parquet')
    report.append(f"\n--- A. {name}: per-gene DEI ({len(piv):,} genes >= {MIN_POS} TpC pos) ---")
    report.append(f"    mean DEI_tpc={piv.tpc.mean():.3e}  mean motif-corrected={piv.mc.mean():.3e}  "
                  f"frac(DEI_tpc>0)={(piv.tpc > 0).mean():.3f}")

from scipy.stats import spearmanr, pearsonr
hA = per_gene_dei([f'BE4_clone{i}' for i in (1, 2, 3, 4)])
hB = per_gene_dei([f'BE4_clone{i}' for i in (5, 6, 7, 8)])
J = hA.join(hB, lsuffix='_A', rsuffix='_B', how='inner')
J = J[(J.npos_tpc_A >= MIN_POS) & (J.npos_tpc_B >= MIN_POS)]
report.append(f"\n--- C. SPLIT-HALF within BE4 (clone1-4 vs clone5-8); {len(J):,} genes ---")
sr = spearmanr(J.tpc_A, J.tpc_B); report.append(f"    DEI_tpc  Spearman={sr.correlation:.3f} (p={sr.pvalue:.1e}) Pearson={pearsonr(J.tpc_A,J.tpc_B)[0]:.3f}")
srm = spearmanr(J.mc_A, J.mc_B); report.append(f"    motif-corrected Spearman={srm.correlation:.3f} (p={srm.pvalue:.1e})")
topA = set(J.sort_values('tpc_A', ascending=False).head(500).index)
topB = set(J.sort_values('tpc_B', ascending=False).head(500).index)
jac = len(topA & topB) / len(topA | topB) if (topA | topB) else 0
report.append(f"    top-500 Jaccard A vs B: {jac:.3f} (random~{500/max(len(J),1):.4f}; enrichment {jac/(500/max(len(J),1)):.1f}x)")
J['qA'] = pd.qcut(J.tpc_A.rank(method='first'), 5, labels=False)
report.append("    mean half-B DEI_tpc by half-A quintile (monotone rise = reproducible):")
report.append("    " + J.groupby('qA').tpc_B.mean().round(7).to_string().replace('\n', '\n    '))

be4 = gene_tables['BE4'].copy()
tot_be4, ed_be4 = pool(COHORTS['BE4']); k = keep & (tot_be4 >= 1)
gsz = pd.Series(genes[k]).groupby(genes[k]).size()
burden = pd.DataFrame({'gene': genes[k], 'ed': ed_be4[k]}).groupby('gene').ed.sum()
be4j = be4.join(gsz.rename('size')).join(burden.rename('burden')).dropna(subset=['size'])
report.append("\n--- E. COVERAGE CONFOUND (BE4) ---")
report.append(f"    DEI_tpc(ratio) vs gene size: Spearman={spearmanr(be4j['size'], be4j['tpc']).correlation:+.3f}  (want ~0)")
report.append(f"    edit-COUNT burden vs gene size: Spearman={spearmanr(be4j['size'], be4j['burden']).correlation:+.3f}  (expect +, confounded)")

b = be4j['burden'].sort_values(ascending=False).values
cum = np.cumsum(b) / b.sum(); ng = len(b)
def frac_genes_for(p): return int(np.searchsorted(cum, p)) / ng
report.append("\n--- F. BURDEN CONCENTRATION (BE4 excess TpC edits) ---")
report.append(f"    50% burden in top {100*frac_genes_for(0.5):.1f}% genes; 80% in {100*frac_genes_for(0.8):.1f}%; 90% in {100*frac_genes_for(0.9):.1f}%")

log('writing per-site label table (kept sites)...')
ks = np.where(keep)[0]
site_df = pd.DataFrame({
    'chrom': chroms[ks], 'pos': positions[ks], 'strand': np.where(strand_plus[ks], '+', '-'),
    'gene': genes[ks], 'motif': np.where(is_tpc[ks], 'TpC', 'nonTpC'),
})
for name, clones in COHORTS.items():
    tot, ed = pool(clones); site_df[f'{name}_ed'] = ed[ks]; site_df[f'{name}_tot'] = tot[ks]
site_df['Parent_ed'] = ped[ks]; site_df['Parent_tot'] = ptot[ks]
site_df.to_parquet(f'{OUT}/per_site_labels_24clone.parquet')
report.append(f"\n--- G. per-site table: {len(site_df):,} sites x {site_df.shape[1]} cols "
              f"(chrom,pos,strand,gene,motif + per-cohort ed/tot) -> per_site_labels_24clone.parquet ---")

txt = '\n'.join(report)
print(txt)
open(f'{OUT}/qc_summary_24clone.txt', 'w').write(txt + '\n')
log('DONE. results in poc24_results/')
