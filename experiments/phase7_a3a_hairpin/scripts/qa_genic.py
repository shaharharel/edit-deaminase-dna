#!/usr/bin/env python
"""Where do the model's top-ranked sites FALL? A coarse version of the spec's tier stratification.

CLAUDE.md's output spec is tier-stratified: Tier A (COSMIC/ClinGen), Tier B (DepMap
essentials), Tier C (coding + UTR + cCRE), Tier D (intergenic). Only refGene is on these
nodes, so Tiers A and B cannot be built tonight -- that needs COSMIC Tier 1, ClinGen HI and
DepMap, and I am not going to approximate a cancer-gene list. What CAN be built is the C/D
split: CDS / UTR / intron / intergenic.

It matters for the safety framing. If the model's top 0.1% is disproportionately INTERGENIC,
a high predicted risk is mostly harmless and the gate would need re-weighting before it means
anything clinically. If it tracks or exceeds the background genic fraction, the ranking is
pointing at places where a mutation can do damage.

COORDINATES, handled explicitly because this is the project's recurring bug: UCSC refGene
gives 0-BASED half-open starts and 1-BASED ends; the trainset's pos is 1-BASED. So a 1-based
position p is inside [txStart, txEnd) iff txStart < p <= txEnd. The script verifies its own
annotation against known genome-wide fractions before reporting anything.
"""
import gzip
import numpy as np
from collections import defaultdict

FEAT = "/data/a3a/feat"
REF = "/data/ref/hg19_refGene.txt.gz"

tx, cds = defaultdict(list), defaultdict(list)
with gzip.open(REF, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c = f[2]
        if not c.startswith("chr") or "_" in c:
            continue
        c = c[3:]
        tx[c].append((int(f[4]), int(f[5])))
        cs, ce = int(f[6]), int(f[7])
        if ce > cs:
            for s, e in zip(f[9].rstrip(",").split(","), f[10].rstrip(",").split(",")):
                s, e = int(s), int(e)
                a, b = max(s, cs), min(e, ce)
                if b > a:
                    cds[c].append((a, b))
print(f"refGene: {sum(len(v) for v in tx.values()):,} transcripts, "
      f"{sum(len(v) for v in cds.values()):,} CDS blocks")

def merge(iv):
    if not iv: return np.zeros((0, 2), np.int64)
    a = np.array(sorted(iv), np.int64)
    out = [a[0].tolist()]
    for s, e in a[1:]:
        if s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return np.array(out, np.int64)

TX = {c: merge(v) for c, v in tx.items()}
CD = {c: merge(v) for c, v in cds.items()}

def inside(iv, pos1):
    """1-based pos is inside 0-based half-open [s,e) iff s < pos <= e."""
    if len(iv) == 0: return np.zeros(len(pos1), bool)
    i = np.searchsorted(iv[:, 0], pos1 - 1, side="right") - 1
    ok = i >= 0
    r = np.zeros(len(pos1), bool)
    r[ok] = pos1[ok] <= iv[i[ok], 1]
    return r

d = np.load(f"{FEAT}/a3a_trainset_v5cov.npz", allow_pickle=True)
y = d["y"].astype(float); chrom = d["chrom"]; pos = d["pos"]
p = np.load(f"{FEAT}/oof_v5cov_hairpin+sequence.npy")
if p.min() < 0 or p.max() > 1: p = 1 / (1 + np.exp(-p))
N = len(y); base = y.mean()

is_tx = np.zeros(N, bool); is_cds = np.zeros(N, bool)
for c in np.unique(chrom):
    m = chrom == c
    is_tx[m] = inside(TX.get(str(c), np.zeros((0, 2), np.int64)), pos[m])
    is_cds[m] = inside(CD.get(str(c), np.zeros((0, 2), np.int64)), pos[m])

cat = np.where(is_cds, "CDS", np.where(is_tx, "intron/UTR", "intergenic"))
print(f"\n=== self-check: are these fractions plausible for hg19? ===")
for k in ("CDS", "intron/UTR", "intergenic"):
    print(f"  background {k:12s} {np.mean(cat == k):.4f}")
print("  expected roughly: CDS ~0.01-0.02, transcribed ~0.40-0.50, intergenic ~0.50-0.60")

k = max(int(round(N * 0.001)), 1)
top = np.argpartition(-p, k - 1)[:k]
print(f"\n=== category composition: background vs the model's top 0.1% (n={k:,}) vs positives ===")
print(f"  {'category':12s} {'background':>11} {'top 0.1%':>10} {'ratio':>7} {'positives':>11} {'ratio':>7}")
for kk in ("CDS", "intron/UTR", "intergenic"):
    b = np.mean(cat == kk)
    t = np.mean(cat[top] == kk)
    q = np.mean(cat[y == 1] == kk)
    print(f"  {kk:12s} {b:>11.4f} {t:>10.4f} {t/b:>6.3f}x {q:>11.4f} {q/b:>6.3f}x")

print(f"\n=== does enrichment differ BY category? (n reported; CDS is thin) ===")
for kk in ("CDS", "intron/UTR", "intergenic"):
    m = cat == kk
    if m.sum() < 5000: 
        print(f"  {kk:12s} n={int(m.sum()):>8,}  (too few)"); continue
    kk2 = max(int(round(m.sum() * 0.001)), 1)
    sub = np.flatnonzero(m)
    t2 = sub[np.argpartition(-p[sub], kk2 - 1)[:kk2]]
    print(f"  {kk:12s} n={int(m.sum()):>8,}  base {y[m].mean():.4f}  "
          f"top0.1% within category {y[t2].mean()/y[m].mean():6.3f}x (n_top={kk2:,})")
