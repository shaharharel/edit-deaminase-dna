#!/usr/bin/env python
"""QA check 4, second half: COMPLEXITY stratification. Never run before.

GC-decile stratification was done (10/10 deciles survived). Sequence COMPLEXITY was
not, and it is the more dangerous of the two here: a hairpin IS an inverted repeat, so
hairpin-rich sites are low-complexity BY CONSTRUCTION. If mutation calls are also
enriched in low-complexity regions - which they are, because alignment is less reliable
there - then "hairpin enrichment" could be a mappability artifact wearing a structural
costume. That is exactly this project's failure family.

Two independent complexity measures on the +-20bp oriented window:
  A. base-composition Shannon entropy (bits, 0-2). Low = AT-rich or homopolymeric.
  B. period-2 self-similarity: fraction of positions with w[i] == w[i+2]. High = a
     dinucleotide repeat like ATATAT or a tandem structure.
For each, decile-stratify and recompute the stem>=6 enrichment WITH an in-stratum
random baseline (200 draws), so the number is never read against a global 1.0.
"""
import numpy as np, sys

P = sys.argv[1] if len(sys.argv) > 1 else '/data/a3a/feat/a3a_trainset_v4.npz'
d = np.load(P, allow_pickle=True)
y, win, stem = d['y'], d['win'], d['stem']
F = win.shape[1] // 2
w = win[:, F - 20:F + 21].astype(np.int8)
base = y.mean()
print(f'=== {P.split("/")[-1]}  rows={len(y):,} pos={int(y.sum()):,} base_rate={base:.4f}')

# A. composition entropy
n = w.shape[1]
p = np.stack([(w == b).sum(axis=1) / n for b in range(4)], axis=1).astype(np.float64)
with np.errstate(divide='ignore', invalid='ignore'):
    ent = -np.nansum(np.where(p > 0, p * np.log2(p), 0.0), axis=1)
# B. period-2 self-similarity
rep2 = (w[:, :-2] == w[:, 2:]).mean(axis=1)

rng = np.random.default_rng(0)
sel = stem >= 6
print(f'  stem>=6 sites: {int(sel.sum()):,}   GLOBAL enrichment '
      f'{y[sel].mean()/base:.4f}x')


def strat(v, label):
    print(f'\n--- {label}')
    print(f'  {"dec":>4} {"range":>13} {"n_tot":>10} {"n_pos_sel":>10} {"enr":>8} {"RANDOM":>8}')
    q = np.quantile(v, np.linspace(0, 1, 11))
    dec = np.clip(np.searchsorted(q[1:-1], v), 0, 9)
    below = 0
    for i in range(10):
        m = dec == i
        s = m & sel
        if s.sum() < 50:
            print(f'  {i:>4} {"":>13} {int(m.sum()):>10,} {"n<50 skipped":>10}')
            continue
        b = y[m].mean()
        enr = y[s].mean() / b
        draws = [y[m][rng.choice(int(m.sum()), int(s.sum()), replace=False)].mean() / b
                 for _ in range(200)]
        r = float(np.mean(draws))
        below += (enr <= r)
        print(f'  {i:>4} {f"{q[i]:.3f}-{q[i+1]:.3f}":>13} {int(m.sum()):>10,} '
              f'{int((y[s]==1).sum()):>10,} {enr:>7.3f}x {r:>7.3f}x')
    print(f'  deciles where stem>=6 FAILS to beat its in-stratum baseline: {below}/10')


strat(ent, 'A. base-composition entropy deciles (low = homopolymeric/AT-rich)')
strat(rep2, 'B. period-2 self-similarity deciles (high = dinucleotide repeat)')

# extreme low-complexity tail: does the signal depend on it?
print('\n--- C. drop the most repetitive sites entirely and re-measure')
for cut in (1.00, 0.95, 0.90, 0.75, 0.50):
    thr = np.quantile(rep2, cut)
    keep = rep2 <= thr
    s = keep & sel
    b = y[keep].mean()
    print(f'  keep bottom {cut*100:5.1f}% by repeat  n={int(keep.sum()):>10,}  '
          f'n_pos_sel={int((y[s]==1).sum()):>7,}  enr={y[s].mean()/b:.4f}x')
