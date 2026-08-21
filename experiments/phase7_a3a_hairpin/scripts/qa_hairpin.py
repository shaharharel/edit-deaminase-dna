#!/usr/bin/env python
"""QA: verify the HAIRPIN FEATURE against hg19 directly.

stem>=6 is what produces the 4.264x tail. Every claim tonight rests on
hairpin_from_windows() having the right orientation and indexing -- and that function
has never been checked against the reference, only trusted. Its self-consistency is
not evidence: bug 2 was an orientation bug that was perfectly self-consistent.

This does NOT reimplement the search (that would just be a second chance to make the
same mistake). It VERIFIES THE CLAIM: for each site the artifact asserts a specific
(stem S, loop L, pos_in_loop k). From the generator:
    ls = C - k ;  le = ls + L
    left  = win[ls-S : ls]
    right = win[le : le+S]
    requires  left == revcomp(right)
    gc_pairs  = #(G or C) in left
    hp_score  = 2*gc + 1*(S-gc) - 0.25*L
So I rebuild the strand-oriented window from hg19 and check those four assertions hold.
A wrong sign, a wrong centre, or a missing reverse-complement all break them.
"""
import numpy as np, sys

REF = '/mnt/data/ref/hg19.fa'
CHROM = sys.argv[1] if len(sys.argv) > 1 else '22'
U = np.load(f'/mnt/data/a3a/feat/universe_chr{CHROM}.npz')
pos, strand = U['pos'], U['strand']
stem, loop, kpos, gcp, hps = U['stem'], U['loop'], U['kpos'], U['gc_pairs'], U['hp_score']

B = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
COMPI = np.array([3, 2, 1, 0, 4], dtype=np.uint8)     # A<->T, C<->G, N->N
seq, on = [], False
for line in open(REF):
    if line[0] == '>':
        n = line[1:].split()[0]
        if on: break
        on = n in (f'chr{CHROM}', CHROM); continue
    if on: seq.append(line.strip())
s = ''.join(seq).upper()
a = np.frombuffer(s.encode(), dtype=np.uint8)
g = np.full(a.shape, 4, dtype=np.uint8)
for ch, v in B.items():
    g[a == ord(ch)] = v
print(f'chr{CHROM} loaded: {len(g):,} bases')

FLANK = 40
rng = np.random.default_rng(0)


def build(idx):
    w = g[np.clip(pos[idx][:, None] + np.arange(-FLANK, FLANK + 1)[None, :], 0, len(g) - 1)]
    neg = strand[idx] == 1
    if neg.any():
        w[neg] = COMPI[w[neg][:, ::-1]]
    return w


def check(idx, label):
    w = build(idx)
    C = FLANK
    S, L, k = stem[idx].astype(int), loop[idx].astype(int), kpos[idx].astype(int)
    ok_pair = np.ones(len(idx), bool)
    ok_gc = np.ones(len(idx), bool)
    ok_sc = np.ones(len(idx), bool)
    ok_focal = w[:, C] == 1                      # focal base must be C after orientation
    for i in range(len(idx)):
        if S[i] == 0:
            continue
        ls = C - k[i]; le = ls + L[i]
        left = w[i, ls - S[i]:ls]; right = w[i, le:le + S[i]]
        ok_pair[i] = bool((left == COMPI[right[::-1]]).all())
        gc = int(((left == 1) | (left == 2)).sum())
        ok_gc[i] = (gc == gcp[idx[i]])
        ok_sc[i] = abs((2.0 * gc + 1.0 * (S[i] - gc) - 0.25 * L[i]) - hps[idx[i]]) < 1e-6
        # focal base must lie inside the loop
        ok_focal[i] &= (ls <= C < le)
    n = len(idx)
    print(f'  {label:22s} n={n:>7,}  stem-pairing {100*ok_pair.mean():7.3f}%  '
          f'gc_pairs {100*ok_gc.mean():7.3f}%  hp_score {100*ok_sc.mean():7.3f}%  '
          f'focal-C-in-loop {100*ok_focal.mean():7.3f}%')
    return ok_pair.all() and ok_gc.all() and ok_sc.all() and ok_focal.all()


print('--- verifying the ASSERTED hairpin against hg19')
allok = True
sel = np.flatnonzero(stem >= 6)
allok &= check(rng.choice(sel, min(20000, len(sel)), replace=False), 'stem>=6 (drives tail)')
sel8 = np.flatnonzero(stem >= 8)
allok &= check(rng.choice(sel8, min(5000, len(sel8)), replace=False), 'stem>=8')
allok &= check(rng.choice(len(pos), 20000, replace=False), 'random sites')
for st in (0, 1):
    m = np.flatnonzero((stem >= 6) & (strand == st))
    allok &= check(rng.choice(m, min(10000, len(m)), replace=False), f'stem>=6 strand={st}')
print(f'\n  VERDICT: {"ALL ASSERTIONS HOLD" if allok else "*** MISMATCH FOUND ***"}')
