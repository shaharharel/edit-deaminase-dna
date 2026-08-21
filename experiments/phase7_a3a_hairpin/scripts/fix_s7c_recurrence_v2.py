#!/usr/bin/env python
"""s7c: exclude editor-specific sites that RECUR across independent clones. (v2)

WHY v2 — a defect in my own v1 patch, caught by reading it before it was ever applied:

v1 used its own marker string as the idempotency check (`if "A3A_EXCLUDE_RECURRENT" in
src`), and the inserted code ALSO contains that string. To stop the check from matching
its own output I obfuscated one occurrence as
    os.environ.get(chr(65)+chr(51)+chr(65)+chr(95)+'EXCLUDE_RECURRENT')
That works, and it is terrible: it writes unreadable gibberish into a script that
someone will read months from now while trying to understand an editor result. A patch
that damages the readability of the thing it patches is not a fix.

v2 uses a distinct sentinel comment for idempotency, so the inserted code can say
exactly what it means.

MEASURED JUSTIFICATION (11:45 UTC): two independent A3A-Y130F clones share 3,000
editor-specific sites against 1.8 expected by chance — 1642x, 15% of each clone's set.
VAF median 0.083, so not lineage germline; at coverage ~30 that is alt=2-3 reads at
recurrent positions. Splitting clone2's sites by whether clone5 also called them:
    stem   shared(n=3,000)   private(n=16,975)
      6        0.795             0.920
      7        1.275             0.994
      8        1.326             1.026
The marginal pooled stem-8 result (1.196, p=0.046) came entirely from the shared
fraction. On clone-private sites the editor is null to three decimals.
"""
import sys, os, shutil

BASE = "/mnt/data/a3a"
p = f"{BASE}/s7c_editor.py"
src = open(p).read()

SENTINEL = "# --- cross-clone recurrence filter (added 2026-08-21) ---"
if SENTINEL in src:
    print("  already patched"); sys.exit(0)

ANCHOR = '    silent = (P["alt"] == 0) & (P["cov"] >= 15)'
if ANCHOR not in src:
    print("  *** ANCHOR NOT FOUND - NOT PATCHED ***"); sys.exit(1)

FILTER = ANCHOR + '''

''' + SENTINEL + '''
    # Two independent A3A-Y130F clones share 3,000 editor-specific sites where chance
    # predicts 1.8 (1642x). Those sites are not clone-private mutation and they carried
    # the whole marginal stem-8 signal (shared 1.326 vs private 1.026). Dropping them is
    # the site-recurrence control that downgraded an earlier phase of this project.
    for _sib in [x for x in RECUR_SIBS if x]:
        _S = np.load(f"{FEAT}/counts_{_sib}_chr{c}.npz")
        _before = int(silent.sum())
        silent = silent & ~((_S["cov"] >= 8) & (_S["alt"] >= 2))
        RECUR_DROPPED[0] += _before - int(silent.sum())'''

src = src.replace(ANCHOR, FILTER, 1)

src = src.replace(
    'burden = {b: [0, 0, 0, 0] for b in BANDS}',
    'RECUR_SIBS = [x for x in os.environ.get("A3A_EXCLUDE_RECURRENT", "").split(",") if x]\n'
    'RECUR_DROPPED = [0]\n'
    'burden = {b: [0, 0, 0, 0] for b in BANDS}', 1)

src = src.replace(
    'print("ENDPOINT A -- BURDEN',
    'if RECUR_SIBS:\n'
    '    print(f"CROSS-CLONE RECURRENCE FILTER ACTIVE: dropped {RECUR_DROPPED[0]:,} sites "\n'
    '          f"also specific in {\',\'.join(RECUR_SIBS)}")\n'
    'else:\n'
    '    print("cross-clone recurrence filter: OFF "\n'
    '          "(set A3A_EXCLUDE_RECURRENT=<sibling,sibling> to enable)")\n'
    'print("ENDPOINT A -- BURDEN', 1)

if not src.startswith("#!") or "import os," not in src:
    src = src.replace("import sys, glob, re", "import os, sys, glob, re", 1)

shutil.copy(p, p + ".prerecur")
open(p, "w").write(src)
compile(open(p).read(), p, "exec")
print("  s7c_editor.py: recurrence filter added, compiles (backup .prerecur)")
print("  default OFF; every run states whether it was active and how many sites it cut")
