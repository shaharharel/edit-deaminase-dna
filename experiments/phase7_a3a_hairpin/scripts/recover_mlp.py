#!/usr/bin/env python
"""Recover the MLP ladder's results from its log, and record honestly why they are log-derived.

The run was OOM-KILLED by the kernel on its final block. "ALL (+seq81)" concatenates
ntv2(1024) + hyena(256) + struct(6) + thermo(10) + seq81(324) = 1,620 dims over 923,989 rows
= 6.0 GB, plus an 80% training copy and a 20% test copy, on top of the source blocks already
resident. On a 29 GB host that is over the line.

I found this exact fragility two ticks ago and fixed it in gpu_fuse2.py -- the CNN ladder --
and did not carry the fix to gpu_fuse.py. That is the third time tonight I have fixed the
copy in front of me instead of the family.

Eleven of twelve blocks completed and printed before the kill. Their numbers are real; only
the kitchen-sink block and the ranked summary were lost. Because json.dump ran only after the
loop, nothing was written -- the same single-save flaw that cost the first encoder run at 96%.
"""
import json, re, sys

log = open("/mnt/a3a/fuse2.log", errors="replace").read()
pat = re.compile(
    r"^\s{2}(\S.*?)\s+AUROC\s+([\d.]+)\+-([\d.]+)\s+d=(\d+)\s+\|\s+"
    r"top0\.1%\s+([\d.]+)x\s+\(rand\s+([\d.]+),\s+n=([\d,]+)\)\s+\|\s+"
    r"top1\.0%\s+([\d.]+)x\s+\(rand\s+([\d.]+),\s+n=[\d,]+\)", re.M)
out = {}
for m in pat.finditer(log):
    out[m.group(1).strip()] = {
        "auroc_mean": float(m.group(2)), "auroc_sd": float(m.group(3)),
        "n_feat": int(m.group(4)),
        "top0.1": float(m.group(5)), "rand0.1": float(m.group(6)),
        "n_top0.1": int(m.group(7).replace(",", "")),
        "top1.0": float(m.group(8)), "rand1.0": float(m.group(9)),
    }
if not out:
    sys.exit("*** parsed nothing -- refusing to write an empty artifact ***")
out["_provenance"] = {
    "source": "recovered from fuse2.log; the run was OOM-killed before json.dump",
    "complete": False,
    "missing_block": "ALL (+seq81) -- 1,620 dims, the block that caused the kill",
    "blocks_recovered": sum(1 for k in out if k != "_provenance"),
    "metric": "pooled out-of-fold top-K% enrichment, 5 folds by held-out CHROMOSOME",
    "gb_baseline_top0.1": 5.280,
}
json.dump(out, open("/mnt/a3a/deaminaformer_ablation.json", "w"), indent=2)
print(f"recovered {len(out)-1} blocks -> /mnt/a3a/deaminaformer_ablation.json")
for k, v in sorted(((k, v) for k, v in out.items() if k != "_provenance"),
                   key=lambda kv: -kv[1]["top0.1"]):
    print(f"  {v['top0.1']:6.3f}x (rand {v['rand0.1']:.3f})  AUROC {v['auroc_mean']:.4f}  "
          f"d={v['n_feat']:<5d} {k}")
