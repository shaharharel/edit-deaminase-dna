#!/usr/bin/env python
"""MY OWN BUG, and it is the standing family: np.savez APPENDS .npz.

My atomic-write patch to s6_pileup.py wrote:
    _tmp = out + ".partial"
    np.savez_compressed(_tmp, ...)     # numpy appends .npz -> writes <out>.partial.npz
    os.replace(_tmp, out)              # renames <out>.partial, which never existed

So every output landed as counts_<sample>_chr<c>.npz.partial.npz. The driver's
`ls counts_*_chr*.npz | wc -l` COUNTED THEM (they do end in .npz) and reported
"clone6 complete", while np.load(counts_<sample>_chr<c>.npz) raised FileNotFoundError
and the calibrator gate reported "counts incomplete". The critical path stopped at
09:28 for a reason that was neither the data nor the calibrator: it was my patch.

A convention correct locally (savez's auto-extension) consumed downstream as if
universal. Same family as the five known bugs.

Blast radius: clone6 only -- the sole sample pileup'd after the patch. Verified by
counting good vs misnamed names per sample.

This script (1) fixes the GENERATOR, (2) VERIFIES each misnamed file actually loads
and has the right length before touching it, (3) renames atomically.
"""
import os, sys, glob, shutil
import numpy as np

BASE = os.environ.get("A3A_BASE", "/mnt/data/a3a")
FEAT = f"{BASE}/feat"

# ---------- 1. FIX THE GENERATOR (not the artifact) ----------
p = f"{BASE}/s6_pileup.py"
src = open(p).read()
OLD = '    _tmp = out + ".partial"'
NEW = ('    # numpy APPENDS .npz when the name does not already end in it, so the temp name\n'
       '    # must carry the extension or savez writes <out>.partial.npz and the os.replace\n'
       '    # below renames a file that does not exist. That silently mis-named every clone6\n'
       '    # chromosome and stopped the critical path at a gate that blamed the data.\n'
       '    _tmp = out + ".partial.npz"')
if OLD in src:
    shutil.copy(p, p + ".prepartialfix")
    src = src.replace(OLD, NEW, 1)
    open(p, "w").write(src)
    compile(open(p).read(), p, "exec")
    print(f"  GENERATOR FIXED: {p} (backup .prepartialfix), compiles")
elif '.partial.npz"' in src:
    print("  generator already fixed")
else:
    print("  *** generator anchor not found -- NOT patched ***"); sys.exit(1)

# ---------- 2. VERIFY, then 3. RENAME ----------
bad = sorted(glob.glob(f"{FEAT}/counts_*_chr*.npz.partial.npz"))
print(f"\n  mis-named files found: {len(bad)}")
if not bad:
    sys.exit(0)

ok, failed = [], []
for f in bad:
    real = f[:-len(".partial.npz")]          # strip BOTH suffixes -> ...chrN.npz
    c = real.split("_chr")[-1].replace(".npz", "")
    try:
        d = np.load(f)
        cov, alt, af, ar = d["cov"], d["alt"], d["alt_fwd"], d["alt_rev"]
        nu = len(np.load(f"{FEAT}/universe_chr{c}.npz")["pos"])
        defects = (len(cov) != nu) + int((alt > cov).sum()) + int((af + ar != alt).sum())
        if defects:
            failed.append((f, f"{defects} integrity defects")); continue
        if os.path.exists(real):
            failed.append((f, "destination already exists -- refusing to overwrite")); continue
        ok.append((f, real, len(cov)))
    except Exception as e:
        failed.append((f, f"{type(e).__name__}: {e}"))

print(f"  verified sound: {len(ok)}   failed verification: {len(failed)}")
for f, why in failed:
    print(f"    REFUSED {os.path.basename(f)}: {why}")
if failed:
    print("\n  *** NOT renaming anything -- fix the failures first ***"); sys.exit(2)

for f, real, n in ok:
    os.replace(f, real)
print(f"  renamed {len(ok)} files; every one loaded with 0 integrity defects first")
print(f"  sites per chromosome checked against universe length in each case")
