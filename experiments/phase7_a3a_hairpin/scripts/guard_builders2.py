#!/usr/bin/env python
"""v2 -- v1 inserted a column-0 guard after a line that was INSIDE an indented
block, producing an IndentationError. compile() caught it before anything ran, and
the .preguard backups make the restore exact. v2 indents the guard to match the
anchor line's own indentation, which is what it should have done to begin with.
"""
import re, shutil, os, textwrap

BODY = '''# --- selection guard (added 2026-08-22 QA) ---
# A threshold correct for one version of an artifact and catastrophic for another
# is this project's recurring bug. MEASURED: tcw>=0.20 & ytca>=0.70 selects 22
# donors under a3a_donor_ranking_v2.tsv and ZERO under v1, whose ytca_frac column
# is corrupt (differs for 1,611/1,636 donors, max |diff| 0.317).
# A collapsed selection must fail loudly, never build an empty trainset.
assert len(sel) >= 10, (
    "selection collapsed to %d donors at TCW_MIN=%s YTCA_MIN=%s -- check WHICH "
    "ranking file was read; v1 gives 0 donors at ytca>=0.70" % (len(sel), TCW_MIN, YTCA_MIN))
print("  selection guard: %d donors at tcw>=%s ytca>=%s" % (len(sel), TCW_MIN, YTCA_MIN))
'''

for f in ("build_A3A.py", "build_A3B.py", "build_A3Awide.py"):
    p = "/data/a3a/" + f
    if not os.path.exists(p):
        print("  %-18s absent" % f); continue
    if os.path.exists(p + ".preguard"):          # undo v1's damage first
        shutil.copy(p + ".preguard", p)
        print("  %-18s restored from .preguard" % f)
    src = open(p).read()
    m = re.search(r"^([ \t]*)sel = \{r\[0\].*$", src, re.M)
    if not m:
        print("  %-18s *** anchor not found, NOT patched ***" % f); continue
    ind = m.group(1)
    guard = "\n" + textwrap.indent(BODY, ind)
    src = src[:m.end()] + guard + src[m.end():]
    shutil.copy(p, p + ".preguard")
    open(p, "w").write(src)
    compile(open(p).read(), p, "exec")           # must pass or we do not keep it
    print("  %-18s guarded at indent %-2d, COMPILES" % (f, len(ind)))
