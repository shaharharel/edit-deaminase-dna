"""Model-guided targeted-panel enrichment demo — f-only (motif) prototype.

Claim: selecting a small region by the model's predicted off-target preference yields a
pooled editing index FAR above the genome-wide index — real (BE4 enriched, nCas9/Parent flat),
out-of-sample (motif spectrum derived from PUBLISHED off-target catalogs, NOT our WGS).

f-only model: rAPOBEC1 trinucleotide (NCN, edited C centered) log-odds, learned from
data/processed/canonical_sites.parquet (CBE/rAPOBEC1 published off-targets) vs the candidate-C
background (our CDS panel). Score every panel C, rank, take top-K%, compute pooled editing index.

Inputs (all local): poc12_results/per_site_labels_12clone.parquet, /tmp/poc_dna/ref/hg19.fa,
data/processed/canonical_sites.parquet. Writes poc12_results/enrichment_demo_fonly.txt.
fxg (accessibility) ablation is a later add once tracks are wired to CDS coords.
"""
import os, sys, time
import numpy as np, pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = os.path.join(REPO, 'experiments/phase4_index/poc12_results')
SITES = f'{RES}/per_site_labels_12clone.parquet'
REF = '/tmp/poc_dna/ref/hg19.fa'
CANON = os.path.join(REPO, 'data/processed/canonical_sites.parquet')
KS = [0.1, 0.5, 1.0, 5.0, 10.0]   # top-K% panels
def log(*a): print(f'[{time.strftime("%H:%M:%S")}]', *a, file=sys.stderr, flush=True)

COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]

# ---------- 1. rAPOBEC1 trinucleotide spectrum from published catalog ----------
log('loading published catalog -> rAPOBEC1 trinuc spectrum...')
cs = pd.read_parquet(CANON)
be4 = cs[(cs.edit_class == 'CBE') & (cs.family == 'rAPOBEC1')]
# normalize trinuc to C-centered uppercase 3mers
tri = be4.trinuc.astype(str).str.upper()
tri = tri[tri.str.len() == 3]
off_counts = tri.value_counts()
log(f'published rAPOBEC1 off-targets: {len(tri):,}; distinct trinuc: {off_counts.size}')

# ---------- 2. load per-site panel ----------
log('loading per-site label table...')
st = pd.read_parquet(SITES, columns=['chrom', 'pos', 'strand', 'gene', 'motif',
                                      'BE4_ed', 'BE4_tot', 'nCas9_ed', 'nCas9_tot',
                                      'Parent_ed', 'Parent_tot'])
N = len(st); log(f'panel positions: {N:,}')

# ---------- 3. trinucleotide context per position from hg19 ----------
log('streaming hg19.fa for trinucleotide context (C-centered, edited strand)...')
by = {}
for i, (c, p, s) in enumerate(zip(st.chrom.values, st.pos.values, st.strand.values)):
    by.setdefault(c, []).append(i)
trinuc = np.empty(N, dtype=object)
want = set(by)
cur = None; parts = None
def emit(chrom, seq):
    L = len(seq)
    for i in by[chrom]:
        p = int(st.pos.values[i]); j = p - 1            # 0-based index of the C/G
        if j - 1 < 0 or j + 1 >= L: trinuc[i] = 'NNN'; continue
        t = seq[j-1:j+2]
        trinuc[i] = t if st.strand.values[i] == '+' else rc(t)
with open(REF) as fh:
    for ln in fh:
        if ln.startswith('>'):
            if cur is not None and parts is not None: emit(cur, ''.join(parts).upper())
            nm = ln[1:].split()[0]; cur = nm; parts = [] if nm in want else None
        elif parts is not None: parts.append(ln.rstrip())
    if cur is not None and parts is not None: emit(cur, ''.join(parts).upper())
st['trinuc'] = trinuc
ok = st.trinuc.str.match('^[ACGT]C[ACGT]$').fillna(False)   # well-formed C-centered
log(f'well-formed C-centered trinuc: {ok.sum():,}/{N:,}')

# ---------- 4. f-only log-odds score ----------
bg_counts = st.loc[ok, 'trinuc'].value_counts()
all_tri = sorted(set(bg_counts.index) | set(off_counts.index))
a = 1.0
off_tot = off_counts.sum() + a * len(all_tri)
bg_tot = bg_counts.sum() + a * len(all_tri)
logodds = {t: np.log((off_counts.get(t, 0) + a) / off_tot) - np.log((bg_counts.get(t, 0) + a) / bg_tot)
           for t in all_tri}
st['f'] = st.trinuc.map(logodds).fillna(-1e9).astype(float)
st.loc[~ok, 'f'] = -1e9                                   # malformed -> never selected
log('top trinuc by f-score: ' + ', '.join(f'{t}:{logodds[t]:+.2f}' for t in sorted(logodds, key=logodds.get, reverse=True)[:6]))

# ---------- 5. editing index over position sets ----------
def EI(mask, ed, tot):
    T = tot[mask].sum(); E = ed[mask].sum()
    return (E / T) if T else 0.0

rng = np.random.default_rng(0)
order = np.argsort(-st.f.values)                          # high f first
cohorts = {'BE4': ('BE4_ed', 'BE4_tot'), 'nCas9': ('nCas9_ed', 'nCas9_tot'), 'Parent': ('Parent_ed', 'Parent_tot')}
ed = {k: st[v[0]].values.astype(np.int64) for k, v in cohorts.items()}
tot = {k: st[v[1]].values.astype(np.int64) for k, v in cohorts.items()}

glob = {k: EI(np.ones(N, bool), ed[k], tot[k]) for k in cohorts}
par_glob = glob['Parent']

rep = []
rep.append("=== MODEL-GUIDED TARGETED-PANEL ENRICHMENT — f-only (rAPOBEC1 trinucleotide) ===")
rep.append(f"panel N={N:,} CDS C/G positions; motif spectrum from {len(tri):,} published rAPOBEC1 off-targets (independent of WGS)")
rep.append(f"\nGENOME-WIDE pooled C->T rate: BE4={glob['BE4']:.3e}  nCas9={glob['nCas9']:.3e}  Parent={glob['Parent']:.3e}")
rep.append(f"GENOME-WIDE BE4 DEI (BE4-Parent) = {glob['BE4']-par_glob:+.3e}")
rep.append("\n--- top-K% model panel vs genome-wide (BE4) and controls ---")
rep.append(f"{'K%':>5} {'#sites':>10} {'#genes':>7} {'BE4 EI':>10} {'BE4 fold':>9} {'BE4 DEI':>10} {'nCas9 fold':>10} {'Parent fold':>11}")
for K in KS:
    n = max(1, int(N * K / 100))
    sel = np.zeros(N, bool); sel[order[:n]] = True
    ngenes = st.loc[sel, 'gene'].nunique()
    be4_top = EI(sel, ed['BE4'], tot['BE4']); par_top = EI(sel, ed['Parent'], tot['Parent']); nc_top = EI(sel, ed['nCas9'], tot['nCas9'])
    be4_fold = be4_top / glob['BE4'] if glob['BE4'] else float('nan')
    nc_fold = nc_top / glob['nCas9'] if glob['nCas9'] else float('nan')
    par_fold = par_top / par_glob if par_glob else float('nan')
    rep.append(f"{K:>5} {n:>10,} {ngenes:>7,} {be4_top:>10.3e} {be4_fold:>9.2f} {be4_top-par_top:>+10.3e} {nc_fold:>10.2f} {par_fold:>11.2f}")

# random-panel control (size-matched), BE4, mean of 5 draws at K=1%
rep.append("\n--- size-matched RANDOM panel control (K=1%), BE4 ---")
n1 = max(1, int(N * 0.01))
rfolds = []
for seed in range(5):
    r = np.random.default_rng(seed); idx = r.choice(N, n1, replace=False)
    m = np.zeros(N, bool); m[idx] = True
    rfolds.append(EI(m, ed['BE4'], tot['BE4']) / glob['BE4'])
rep.append(f"  random 1% BE4 fold = {np.mean(rfolds):.2f} +/- {np.std(rfolds):.2f}  (vs model top-1% above)")

# TpC-only baseline (what naive motif gating gives), for context
tpc = (st.motif.values == 'TpC')
rep.append(f"\n--- naive TpC-gate baseline (no ranking): {tpc.sum():,} sites ({100*tpc.mean():.1f}%) ---")
rep.append(f"  BE4 EI(TpC)={EI(tpc, ed['BE4'], tot['BE4']):.3e}  fold={EI(tpc, ed['BE4'], tot['BE4'])/glob['BE4']:.2f}  (model concentrates beyond flat TpC gate)")

txt = '\n'.join(rep)
print(txt)
open(f'{RES}/enrichment_demo_fonly.txt', 'w').write(txt + '\n')
st[['chrom', 'pos', 'strand', 'gene', 'motif', 'trinuc', 'f']].to_parquet(f'{RES}/site_fscore_fonly.parquet')
log('DONE -> enrichment_demo_fonly.txt + site_fscore_fonly.parquet')
