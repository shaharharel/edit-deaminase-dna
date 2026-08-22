#!/usr/bin/env python
"""LEAVE-ONE-DONOR-OUT error bar on the 5.280x that carries this project's headline.

WHY. v5 is the trainset behind the transfer test's positive control, the architecture
recommendation, and the 5.280x baseline -- and it rests on 21 DONORS. Every enrichment I
have quoted from it carries a chromosome-fold error bar, which measures spatial stability,
NOT donor stability. If one donor's mutations carry the tail, 5.280x is a fact about one
person. Nothing published so far has tested that.

DESIGN. Jackknife: drop donor d entirely, re-run the SAME 5-fold-held-out-CHROMOSOME OOF,
recompute top-0.1% enrichment. 21 replicates plus the full-data reference. Feature block is
hairpin+sequence (86), the exact block that produced 5.280x.

The spread across replicates is the honest error bar. Compare it against the project's
measured ~12% random-baseline noise floor: if the jackknife spread is inside the floor, the
number is donor-robust; if one donor's removal moves it more than the floor, say so.
"""
import numpy as np, time, sys
from sklearn.ensemble import HistGradientBoostingClassifier

FEAT, TS, SEED, CTX = "/data/a3a/feat", "a3a_trainset_v5.npz", 0, 10
K = 0.1  # top-K% -- the tail is the endpoint, not AUROC
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

d = np.load(f"{FEAT}/{TS}", allow_pickle=True)
y = d["y"].astype(np.int8); chrom = d["chrom"]; donor = d["donor"]
hp = np.column_stack([d["stem"], d["loop"], d["kpos"], d["gc_pairs"],
                      d["hp_score"], d["local_gc"]]).astype(np.float32)
win = d["win"]; F = win.shape[1] // 2
offs = [o for o in range(-CTX, CTX + 1) if o != 0]
sq = np.zeros((len(y), len(offs) * 4), dtype=np.float32)
for j, o in enumerate(offs):
    col = win[:, F + o]
    for b in range(4):
        sq[:, j * 4 + b] = (col == b)
X = np.hstack([hp, sq])
log(f"X {X.shape}  positives {int(y.sum()):,}  donors {len(np.unique(donor)):,}")

chroms = sorted(set(chrom.tolist()), key=lambda c: (len(c), c))
folds = [chroms[i::5] for i in range(5)]
rng = np.random.default_rng(SEED)

def oof_enrichment(mask):
    """5-fold OOF by held-out chromosome on the subset, then top-K% enrichment."""
    Xs, ys, cs = X[mask], y[mask], chrom[mask]
    oof = np.zeros(len(ys), dtype=np.float32)
    for hold in folds:
        te = np.isin(cs, hold); tr = ~te
        if te.sum() == 0 or ys[tr].sum() == 0:
            continue
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1,
                early_stopping=True, validation_fraction=0.1, random_state=SEED)
        m.fit(Xs[tr], ys[tr])
        oof[te] = m.predict_proba(Xs[te])[:, 1]
    n = max(1, int(len(ys) * K / 100))
    base = ys.mean()
    top = np.argsort(-oof)[:n]
    enr = ys[top].mean() / base
    sh = rng.permutation(len(ys))[:n]          # RANDOM baseline, same n
    rnd = ys[sh].mean() / base
    return enr, rnd, n, int(ys.sum())

log("=== full data reference (must land near 5.280x) ===")
e, r, n, npos = oof_enrichment(np.ones(len(y), bool))
log(f"  FULL   enr={e:6.3f}x  RANDOM={r:5.3f}x  n={n:,}  pos={npos:,}")
full = e

donors = np.unique(donor)
res = []
for i, dn in enumerate(donors):
    m = donor != dn
    e, r, n, npos = oof_enrichment(m)
    dropped = int((donor == dn).sum())
    res.append((str(dn), e, r, dropped))
    log(f"  drop {str(dn)[:18]:18s} enr={e:6.3f}x  RANDOM={r:5.3f}x  (dropped {dropped:,} rows)")

vals = np.array([x[1] for x in res])
print("\n" + "=" * 66)
print(f"  full-data value                {full:.3f}x")
print(f"  jackknife mean                 {vals.mean():.3f}x")
print(f"  jackknife sd                   {vals.std(ddof=1):.3f}")
print(f"  jackknife min / max            {vals.min():.3f} / {vals.max():.3f}")
print(f"  spread (max-min)               {vals.ptp():.3f}")
n_d = len(vals)
se = np.sqrt((n_d - 1) / n_d * ((vals - vals.mean()) ** 2).sum())
print(f"  jackknife SE on the estimate   {se:.3f}")
print(f"\n  project random-baseline noise floor is ~12% of the value = {0.12*full:.3f}")
print(f"  spread / floor = {vals.ptp()/(0.12*full):.2f}x   ->  "
      f"{'DONOR-ROBUST' if vals.ptp() < 0.12*full else 'DONOR-SENSITIVE, say so'}")
worst = res[int(np.argmin(vals))]
print(f"  most influential donor: {worst[0]} -- removing it gives {worst[1]:.3f}x "
      f"({worst[1]-full:+.3f} vs full)")
