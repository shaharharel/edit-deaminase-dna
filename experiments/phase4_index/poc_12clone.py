"""12-clone cross-deaminase DNA editing-index POC (Doman 2020, hg19).

Larger successor to poc_dna_editing_index.py / poc_splithalf.py. Now genuinely
cross-deaminase: pools COHORTS (BE4 vs YE1-BE4 vs nCas9 vs Parent) rather than a
single editor, so we can test the discriminators that a single clone cannot
(the nCas9 clonal-somatic confound: see memory ncas9-clonal-somatic-confound).

Cohorts (12 clones in home as of build time):
  BE4      = clone1,2,3,5,6,7  (rAPOBEC1, dirty)        n=6
  YE1-BE4  = clone4,5,7        (engineered safer)        n=3
  nCas9    = clone1,6          (deaminase-free control)  n=2
  Parent   = Parent_WGS        (untreated baseline)      n=1

Levanon-faithful: pooled (read-summed) RATIO estimator, NO per-site coverage floor,
germline whitelist via Parent VAF>0.05. Motif channels TpC (APOBEC1 sig) vs non-TpC
(per-gene internal negative control). DEI = treated_ratio - parent_ratio.

Analyses:
  A. Per-cohort pooled per-gene DEI (tpc, npc, motif-corrected=tpc-npc).
  B. Cross-deaminase BULK ordering (pooled genome-wide TpC C->T rate per cohort).
  C. Split-half reproducibility WITHIN BE4 (3v3: {1,2,3} vs {5,6,7}) per-gene DEI Spearman.
  D. Cohort-vs-cohort MOTIF SPECIFICITY (the key discriminator): pooled (TpC-rate / nonTpC-rate)
     per cohort. Hypothesis: BE4 motif-specificity >> nCas9 even though nCas9 bulk C->T is
     confounded by endogenous-APOBEC somatic background.
  E. Coverage confound check: per-gene DEI vs gene size (expect ~0); count-burden vs size (expect +).
  F. Burden concentration / precision@K per cohort.
  G. Per-site label table retaining chrom,pos,strand,gene,motif (+ per-cohort edited/total)
     so Evo context windows can be regenerated later for the Phase-2 model.

Reads /tmp/poc_dna/{mpileups,bed}. Writes experiments/phase4_index/poc12_results/.
"""
import os, sys, time
import numpy as np, pandas as pd

DATA = '/tmp/poc_dna'; MP = f'{DATA}/mpileups'; BED = f'{DATA}/bed'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'poc12_results')
os.makedirs(OUT, exist_ok=True)
os.makedirs(f'{OUT}/figures', exist_ok=True)

COHORTS = {
    'BE4':     ['BE4_clone1', 'BE4_clone2', 'BE4_clone3', 'BE4_clone5', 'BE4_clone6', 'BE4_clone7'],
    'YE1-BE4': ['YE1-BE4_clone4', 'YE1-BE4_clone5', 'YE1-BE4_clone7'],
    'nCas9':   ['nCas9_clone1', 'nCas9_clone6'],
}
PARENT = 'Parent_WGS'
ALL_CLONES = [c for v in COHORTS.values() for c in v] + [PARENT]
PARENT_VAF_CAP = 0.05
MIN_POS = 20   # genes need >=20 TpC positions for per-gene stats

def log(*a): print(f'[{time.strftime("%H:%M:%S")}]', *a, file=sys.stderr, flush=True)

# ---------------- load BED -> position index ----------------
log('loading BED...')
pos_index = {}; chroms = []; positions = []; genes = []; motifs = []; strands = []
for bedp, mot in [(f'{BED}/cds_tpc.bed', 'tpc'), (f'{BED}/cds_npc.bed', 'npc')]:
    with open(bedp) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6: continue
            chrom = f[0]; pos = int(f[1]) + 1  # BED 0-based start -> 1-based mpileup pos
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
    """Return (tot, ed) int32 arrays over N positions. edit = C->T equivalent:
    +strand refC -> count T ; -strand refG -> count A."""
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

log('parsing all clones (this is the slow part)...')
CACHE = f'{DATA}/parse_cache_12clone.npz'
S = {}
if os.path.exists(CACHE):
    log(f'loading parse cache {CACHE} ...')
    z = np.load(CACHE, allow_pickle=False)
    if int(z['N']) == N and all(f'{cl}__tot' in z.files for cl in ALL_CLONES):
        for cl in ALL_CLONES:
            S[cl] = (z[f'{cl}__tot'], z[f'{cl}__ed'])
        log('parse cache OK -> skipped 22min parse')
    else:
        log('parse cache stale (N or clone set changed); reparsing')
        S = {}
if not S:
    for cl in ALL_CLONES:
        S[cl] = parse_sample(cl)
    np.savez(CACHE, N=np.int64(N),
             **{f'{cl}__tot': S[cl][0] for cl in ALL_CLONES},
             **{f'{cl}__ed': S[cl][1] for cl in ALL_CLONES})
    log(f'saved parse cache -> {CACHE}')
ptot, ped = S[PARENT]
parent_vaf = np.divide(ped, ptot, out=np.zeros(N), where=ptot > 0)
keep = (ptot >= 1) & (parent_vaf <= PARENT_VAF_CAP)   # germline whitelist
log(f'kept after germline filter: {keep.sum():,}/{N:,}')

def pool(clones):
    tot = np.zeros(N, dtype=np.int64); ed = np.zeros(N, dtype=np.int64)
    for cl in clones: tot += S[cl][0]; ed += S[cl][1]
    return tot, ed

# ---------------- B. cross-deaminase BULK ordering ----------------
def bulk_rate(clones, mask):
    tot, ed = pool(clones); m = keep & mask
    T = tot[m].sum(); E = ed[m].sum()
    return (E / T if T else 0.0), int(E), int(T)

report = []
report.append("=== 12-CLONE CROSS-DEAMINASE DNA EDITING-INDEX POC ===")
report.append(f"positions N={N:,}; kept(germline)={int(keep.sum()):,}")
report.append("\n--- B. BULK pooled C->T rate (TpC channel), per cohort ---")
prate, pe, pt = bulk_rate([PARENT], is_tpc)
report.append(f"  Parent : {prate:.3e}  ({pe:,}/{pt:,})")
for name, clones in COHORTS.items():
    r, e, t = bulk_rate(clones, is_tpc)
    report.append(f"  {name:8s}: {r:.3e}  ({e:,}/{t:,})  vs Parent x{r/prate:.3f}")

# ---------------- D. cohort-vs-cohort MOTIF SPECIFICITY ----------------
report.append("\n--- D. MOTIF SPECIFICITY: pooled TpC rate / nonTpC rate (>1 = APOBEC-like) ---")
for name, clones in ([('Parent', [PARENT])] + list(COHORTS.items())):
    rt, _, _ = bulk_rate(clones, is_tpc)
    rn, _, _ = bulk_rate(clones, ~is_tpc)
    spec = rt / rn if rn else float('nan')
    report.append(f"  {name:8s}: TpC={rt:.3e} nonTpC={rn:.3e}  specificity={spec:.3f}  (TpC-nonTpC={rt-rn:+.3e})")

# ---------------- A. per-cohort per-gene DEI ----------------
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
    npos = r[r.tag == 'tpc'].set_index('gene').npos
    piv['npos_tpc'] = npos
    piv['mc'] = piv['tpc'] - piv['npc']  # motif-corrected
    return piv

log('per-gene DEI per cohort...')
gene_tables = {}
for name, clones in COHORTS.items():
    piv = per_gene_dei(clones)
    piv = piv[piv.npos_tpc >= MIN_POS]
    gene_tables[name] = piv
    piv.to_parquet(f'{OUT}/per_gene_dei_{name}.parquet')
    report.append(f"\n--- A. {name}: per-gene DEI ({len(piv):,} genes >= {MIN_POS} TpC pos) ---")
    report.append(f"    mean DEI_tpc={piv.tpc.mean():.3e}  mean motif-corrected={piv.mc.mean():.3e}  "
                  f"frac(DEI_tpc>0)={(piv.tpc > 0).mean():.3f}")

# ---------------- C. split-half within BE4 (3 v 3) ----------------
from scipy.stats import spearmanr, pearsonr
hA = per_gene_dei(['BE4_clone1', 'BE4_clone2', 'BE4_clone3'])
hB = per_gene_dei(['BE4_clone5', 'BE4_clone6', 'BE4_clone7'])
J = hA.join(hB, lsuffix='_A', rsuffix='_B', how='inner')
J = J[(J.npos_tpc_A >= MIN_POS) & (J.npos_tpc_B >= MIN_POS)]
report.append(f"\n--- C. SPLIT-HALF within BE4 (clone1,2,3 vs clone5,6,7); {len(J):,} genes ---")
sr = spearmanr(J.tpc_A, J.tpc_B); report.append(f"    DEI_tpc  Spearman={sr.correlation:.3f} (p={sr.pvalue:.1e}) Pearson={pearsonr(J.tpc_A,J.tpc_B)[0]:.3f}")
srm = spearmanr(J.mc_A, J.mc_B); report.append(f"    motif-corrected Spearman={srm.correlation:.3f} (p={srm.pvalue:.1e})")
topA = set(J.sort_values('tpc_A', ascending=False).head(500).index)
topB = set(J.sort_values('tpc_B', ascending=False).head(500).index)
jac = len(topA & topB) / len(topA | topB) if (topA | topB) else 0
report.append(f"    top-500 Jaccard A vs B: {jac:.3f} (random~{500/max(len(J),1):.4f}; enrichment {jac/(500/max(len(J),1)):.1f}x)")
J['qA'] = pd.qcut(J.tpc_A.rank(method='first'), 5, labels=False)
report.append("    mean half-B DEI_tpc by half-A quintile (monotone rise = reproducible):")
report.append("    " + J.groupby('qA').tpc_B.mean().round(7).to_string().replace('\n', '\n    '))

# ---------------- E. coverage confound check (BE4) ----------------
be4 = gene_tables['BE4'].copy()
# gene size proxy = TpC positions; ratio DEI should be ~0 vs size, count-burden should be +
tot_be4, ed_be4 = pool(COHORTS['BE4'])
k = keep & (tot_be4 >= 1)
gsz = pd.Series(genes[k]).groupby(genes[k]).size()
burden = pd.DataFrame({'gene': genes[k], 'ed': ed_be4[k]}).groupby('gene').ed.sum()
be4j = be4.join(gsz.rename('size')).join(burden.rename('burden')).dropna(subset=['size'])
csz_dei = spearmanr(be4j['size'], be4j['tpc'])
csz_burden = spearmanr(be4j['size'], be4j['burden'])
report.append("\n--- E. COVERAGE CONFOUND (BE4) ---")
report.append(f"    DEI_tpc(ratio) vs gene size: Spearman={csz_dei.correlation:+.3f}  (want ~0)")
report.append(f"    edit-COUNT burden vs gene size: Spearman={csz_burden.correlation:+.3f}  (expect +, confounded)")

# ---------------- F. burden concentration (BE4) ----------------
b = be4j['burden'].sort_values(ascending=False).values
cum = np.cumsum(b) / b.sum()
ng = len(b)
def frac_genes_for(p): return int(np.searchsorted(cum, p)) / ng
report.append("\n--- F. BURDEN CONCENTRATION (BE4 excess TpC edits) ---")
report.append(f"    50% burden in top {100*frac_genes_for(0.5):.1f}% genes; 80% in {100*frac_genes_for(0.8):.1f}%; 90% in {100*frac_genes_for(0.9):.1f}%")

# ---------------- G. per-site label table (coords retained for Evo) ----------------
log('writing per-site label table (kept sites)...')
ks = np.where(keep)[0]
site_df = pd.DataFrame({
    'chrom': chroms[ks], 'pos': positions[ks], 'strand': np.where(strand_plus[ks], '+', '-'),
    'gene': genes[ks], 'motif': np.where(is_tpc[ks], 'TpC', 'nonTpC'),
})
for name, clones in COHORTS.items():
    tot, ed = pool(clones)
    site_df[f'{name}_ed'] = ed[ks]; site_df[f'{name}_tot'] = tot[ks]
site_df['Parent_ed'] = ped[ks]; site_df['Parent_tot'] = ptot[ks]
site_df.to_parquet(f'{OUT}/per_site_labels_12clone.parquet')
report.append(f"\n--- G. per-site table: {len(site_df):,} sites x {site_df.shape[1]} cols "
              f"(chrom,pos,strand,gene,motif + per-cohort ed/tot) -> per_site_labels_12clone.parquet ---")

txt = '\n'.join(report)
print(txt)
open(f'{OUT}/qc_summary_12clone.txt', 'w').write(txt + '\n')
log('DONE. results in poc12_results/')
