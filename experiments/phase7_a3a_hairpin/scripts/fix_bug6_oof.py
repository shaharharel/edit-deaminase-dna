#!/usr/bin/env python
"""BUG 6, third sweep: the .npy out-of-fold arrays were never covered by the first two.

My earlier sweeps searched for *.json outputs and declared the family closed. Six scripts
write OOF arrays too, and every one of them hardcodes the v2 name:

    s3_train_v2.py  IN v2      -> oof_v2_{block}.npy
    s3_v4.py        IN v4      -> oof_v2_{block}.npy
    s3_v4gc.py      IN v4      -> oof_v2_{block}.npy
    s3_v4s.py       IN v4s     -> oof_v2_{block}.npy
    s3_v5.py        IN v5      -> oof_v2_{block}.npy
    s3_v5cov.py     IN v5cov   -> oof_v2_{block}.npy

DAMAGE, established by file size rather than by assumption: the surviving files are
7,160,560 bytes = 895,070 float64, which is v5cov's row count (895,054), and their mtimes
19:00-19:18 match the v5cov run. So each run overwrote the last and only v5cov's arrays
remain, filed under 'v2'. The published findings are unaffected -- enrichment and AUROC live
in the correctly named JSONs -- but every earlier run's OOF predictions are gone, and anyone
recomputing a metric from oof_v2_*.npy would have silently used v5cov.

Lesson recorded: 'I fixed that bug family' is a claim about a SEARCH, and the search was
incomplete twice. This sweep covers every output extension, not just the one that bit me.

Fixes the generator in all six, and renames the surviving artifacts to what they actually are.
"""
import sys, os, re, glob, shutil

BASE = sys.argv[1] if len(sys.argv) > 1 else "/data/a3a"
FEAT = f"{BASE}/feat"

for name in sorted(os.path.basename(p) for p in glob.glob(f"{BASE}/s[35]_*.py")):
    p = os.path.join(BASE, name)
    src = open(p).read()
    if 'oof_v2_' not in src:
        continue
    m = re.search(r'(a3a_trainset_[a-z0-9]+\.npz)', src)
    if not m:
        print(f"  {name}: writes oof_v2_ but reads no trainset - NOT PATCHED"); continue
    tag = m.group(1).replace("a3a_trainset_", "").replace(".npz", "")
    if "TAG = TRAINSET" not in src:
        print(f"  {name}: has no TAG constant - NOT PATCHED (fix the json sweep first)"); continue
    new = src.replace('f"{FEAT}/oof_v2_{bname}.npy"', 'f"{FEAT}/oof_{TAG}_{bname}.npy"')
    if new == src:
        print(f"  {name}: oof_v2_ present but not in the expected f-string - NOT PATCHED"); continue
    shutil.copy(p, p + ".preoof")
    open(p, "w").write(new)
    compile(open(p).read(), p, "exec")
    print(f"  {name}: reads {tag} -> now writes oof_{tag}_<block>.npy (backup .preoof)")

# rename the surviving artifacts to what they demonstrably are
import numpy as np
sizes = {}
for f in glob.glob(f"{FEAT}/a3a_trainset_v*.npz"):
    t = os.path.basename(f).replace("a3a_trainset_", "").replace(".npz", "")
    sizes[len(np.load(f, allow_pickle=True)["y"])] = t
print("\n  identifying the orphaned oof_v2_*.npy by row count:")
for f in sorted(glob.glob(f"{FEAT}/oof_v2_*.npy")):
    n = len(np.load(f, mmap_mode="r"))
    owner = sizes.get(n)
    if owner is None:
        print(f"    {os.path.basename(f)}: n={n:,} matches NO trainset - left alone"); continue
    dst = f.replace("oof_v2_", f"oof_{owner}_")
    if os.path.exists(dst):
        print(f"    {os.path.basename(f)}: n={n:,} -> {owner}, but target exists - left alone"); continue
    os.rename(f, dst)
    print(f"    {os.path.basename(f)}: n={n:,} -> {os.path.basename(dst)}")
