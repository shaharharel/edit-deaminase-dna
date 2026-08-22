#!/usr/bin/env python
"""QA of MY OWN calibration claim from last tick.

I reported ECE 0.0013 against a spec of 0.05 and called it "met by a factor of 10-40". Before
that goes any further it needs the control every calibration number needs and that I did not
run: WHAT IS THE ECE OF A MODEL WITH NO SKILL?

A constant predictor that outputs the base rate for every site is perfectly calibrated by
construction -- ECE ~ 0 -- and useless. Since these models barely spread their predictions
(most mass between 0.04 and 0.19 around a base rate of 0.094), a large part of that 0.0013
may be inherited from predicting near-constant rather than earned. Low ECE is necessary, not
sufficient, and reporting it without this control overstates the result.

Also measures calibration AT THE OPERATING POINT, which the quantile-bin curve hides: the top
bin held 6.7% of the data, so the top 0.1% where a gate would sit was diluted inside it.
"""
import numpy as np

FEAT = "/data/a3a/feat"
y = np.load(f"{FEAT}/a3a_trainset_v5cov.npz", allow_pickle=True)["y"].astype(float)
p = np.load(f"{FEAT}/oof_v5cov_hairpin+sequence.npy")
if p.min() < 0 or p.max() > 1:
    p = 1 / (1 + np.exp(-p))
N = len(y); base = y.mean()
print(f"n={N:,}  base rate {base:.5f}")

def ece(pp, yy, nb=15):
    e = np.unique(np.quantile(pp, np.linspace(0, 1, nb + 1)))
    tot = 0.0
    for i in range(len(e) - 1):
        m = (pp >= e[i]) & (pp < e[i + 1]) if i < len(e) - 2 else (pp >= e[i])
        if m.sum() < 50: continue
        tot += m.sum() / len(pp) * abs(pp[m].mean() - yy[m].mean())
    return tot

rng = np.random.default_rng(0)
print(f"\n=== the control I skipped: what does NO SKILL score? ===")
CONTROLS = {
    "constant = base rate (zero skill)": np.full(N, base),
    "base rate + tiny noise":            base + rng.normal(0, 1e-3, N),
    "random uniform in the model range": rng.uniform(p.min(), p.max(), N),
    "the model itself":                  p,
}
for name, pp in CONTROLS.items():
    pp = np.clip(pp, 1e-6, 1 - 1e-6)
    print(f"  {name:36s} Brier {np.mean((pp-y)**2):.5f}   ECE {ece(pp, y):.5f}")
print("  A constant predictor is PERFECTLY calibrated and completely useless. ECE alone")
print("  cannot distinguish it from a good model; that is why it needs a skill metric beside")
print("  it, and why 'ECE 10-40x under spec' was an overstatement on its own.")

print(f"\n=== spread: how much does the model actually move? ===")
q = np.quantile(p, [0.001, 0.01, 0.5, 0.99, 0.999])
print(f"  prediction quantiles 0.1/1/50/99/99.9%: " + " ".join(f"{v:.4f}" for v in q))
print(f"  ratio p99.9 / p0.1 = {q[-1]/q[0]:.2f}x   (a constant predictor would be 1.00x)")

print(f"\n=== calibration AT THE OPERATING POINT (what the quantile curve hid) ===")
for pct in (0.1, 0.5, 1.0, 5.0):
    k = max(int(round(N * pct / 100)), 1)
    idx = np.argpartition(-p, k - 1)[:k]
    pred, obs = p[idx].mean(), y[idx].mean()
    print(f"  top {pct:>4}%  n={k:>7,}  predicted {pred:.4f}  observed {obs:.4f}  "
          f"gap {obs-pred:+.4f}  enrichment {obs/base:5.3f}x")
print("  The gate operates in these rows, not in a 60,000-site quantile bin.")
