#!/usr/bin/env python
"""QA: editor-specific site set — VAF profile and cross-clone recurrence.  (v2, cleaned)

v1 of this script produced correct numbers but carried three defects that would bite the
next person to run it. It is saved on the node, so it is worth fixing rather than
leaving as-is:

  1. DEAD CODE: `spec_sets[s].append(np.flatnonzero(spec) + hash(c) % 1)`.
     `hash(c) % 1` is ALWAYS 0, so the term did nothing; `spec_sets` was never read; and
     it accumulated ~40k-element arrays per chromosome for no purpose.
  2. FRAGILE INIT: `overlap` was created inside the loop under `if c == CH[0]`. That
     works only if CH[0] is genuinely processed first. Reorder CH, or run a subset, and
     it raises NameError mid-run — after minutes of work.
  3. `overlap.setdefault("both", [])` papered over (2) instead of initialising cleanly.

None of this changed the reported result: the overlap arithmetic never touched
`spec_sets`, and CH[0] was processed first. But an instrument that only works when its
inputs arrive in one particular order is exactly the class of thing that has misled me
seven times tonight, so it gets fixed rather than trusted.
"""
import numpy as np

FEAT = "/mnt/data/a3a/feat"
CH = [str(i) for i in range(1, 23)] + ["X"]
ED = ["P66-A3A-Y130F-clone2", "P66-A3A-Y130F-clone5"]
CAL = "P66-D10A-clone1"
BINS = [(0.0, 0.05), (0.05, 0.15), (0.15, 0.35), (0.35, 1.01)]

prof = {s: np.zeros(len(BINS), dtype=np.int64) for s in ED + [CAL]}
nspec = {s: 0 for s in ED + [CAL]}
ov = {"a": 0, "b": 0, "both": 0, "elig": 0}          # initialised ONCE, before the loop

for c in CH:
    P = np.load(f"{FEAT}/counts_Parent_chr{c}.npz")
    silent = (P["alt"] == 0) & (P["cov"] >= 15)

    for s in ED + [CAL]:
        D = np.load(f"{FEAT}/counts_{s}_chr{c}.npz")
        cov, alt = D["cov"], D["alt"]
        spec = silent & (cov >= 8) & (alt >= 2)
        nspec[s] += int(spec.sum())
        v = np.zeros(len(cov))
        nz = cov > 0
        v[nz] = alt[nz] / cov[nz]
        vv = v[spec]
        for i, (lo, hi) in enumerate(BINS):
            prof[s][i] += int(((vv >= lo) & (vv < hi)).sum())

    A = np.load(f"{FEAT}/counts_{ED[0]}_chr{c}.npz")
    B = np.load(f"{FEAT}/counts_{ED[1]}_chr{c}.npz")
    el = silent & (A["cov"] >= 8) & (B["cov"] >= 8)   # JOINTLY eligible: the right base
    sa = el & (A["alt"] >= 2)
    sb = el & (B["alt"] >= 2)
    ov["elig"] += int(el.sum())
    ov["a"] += int(sa.sum())
    ov["b"] += int(sb.sum())
    ov["both"] += int((sa & sb).sum())

print("=== A. VAF profile of the analysis set (Parent-masked, cov>=8, alt>=2) ===")
print(f"  {'sample':24s} {'n_spec':>10s} {'<0.05':>8s} {'.05-.15':>8s} {'.15-.35':>8s} {'>=0.35':>8s}")
for s in ED + [CAL]:
    t = max(prof[s].sum(), 1)
    print(f"  {s:24s} {nspec[s]:>10,} " + " ".join(f"{100 * x / t:7.1f}%" for x in prof[s]))

print("\n=== B. cross-clone recurrence of editor-specific sites ===")
exp = ov["a"] * ov["b"] / max(ov["elig"], 1)
print(f"  jointly eligible   {ov['elig']:>12,}")
print(f"  clone2 specific    {ov['a']:>12,}")
print(f"  clone5 specific    {ov['b']:>12,}")
print(f"  shared by BOTH     {ov['both']:>12,}")
print(f"  expected by chance {exp:>12,.1f}")
print(f"  OBS / EXP          {ov['both'] / max(exp, 1e-9):>12.2f}x")
print("  ~1x => clone-private events.  >>1x => a systematic artefact shared by both,")
print("  which must be excluded before Endpoint B is computed.")
