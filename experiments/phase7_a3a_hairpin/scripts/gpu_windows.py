#!/usr/bin/env python
"""Extract EXTENDED-CONTEXT windows for the GPU architecture track.

The GB baseline sees 81 bp. Its ceiling is what 81 bp can carry, and the structure branch
already showed that information ABOUT the window (hairpin geometry) beats more of the window
itself. The architecture question is whether long-range context -- 1 kb, where replication
timing, nucleosome phasing and larger secondary structure live -- adds anything the 81 bp
window cannot express. A DNA language model is the natural encoder for that; this script
just produces its input.

+-512 bp around each site, strand-oriented exactly as the 81 bp windows were, written as
uint8 codes A/C/G/T/N = 0/1/2/3/4.

VERIFICATION IS BUILT IN, not optional: the extracted 1 kb window must agree with the
existing 81 bp window over their shared span, base for base. A silent coordinate or strand
error here would be bug 3 all over again, and it would be invisible in the trained model --
it would just look like the architecture failing.
"""
import sys, os
import numpy as np

BASE = "/data/a3a"
FEAT = f"{BASE}/feat"
SRC = sys.argv[1] if len(sys.argv) > 1 else "a3a_trainset_v5.npz"
HALF = 512
OUT = f"{FEAT}/{SRC.replace('.npz','')}_win1k.npz"

import pyfaidx
REF = None
for p in ("/data/ref/hg19.fa", f"{BASE}/ref/hg19.fa", "/mnt/data/a3a/ref/hg19.fa"):
    if os.path.exists(p):
        REF = p; break
if REF is None:
    raise SystemExit("FATAL: hg19.fa not found in /data/ref, {BASE}/ref, /mnt/data/a3a/ref")
print(f"reference: {REF}")
# /data/ref is SHARED with the RNA sister project. Keep the .fai index inside /data/a3a so
# nothing is written outside our own tree, even a harmless index file.
IDX = f"{FEAT}/hg19.fa.fai.a3a"
fa = pyfaidx.Fasta(REF, as_raw=True, sequence_always_upper=True, indexname=IDX)

d = np.load(f"{FEAT}/{SRC}", allow_pickle=True)
chrom, pos, strand, win81 = d["chrom"], d["pos"], d["strand"], d["win"]
N = len(pos)
print(f"{SRC}: n={N:,}  extracting +-{HALF} bp")

CODE = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate("ACGT"):
    CODE[ord(b)] = i
COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)

out = np.full((N, 2 * HALF + 1), 4, dtype=np.uint8)
keys = {k.split()[0]: k for k in fa.keys()}
mid81 = win81.shape[1] // 2
for i in range(N):
    c = str(chrom[i])
    k = keys.get(c) or keys.get("chr" + c)
    if k is None:
        continue
    p0 = int(pos[i]) - 1
    lo, hi = p0 - HALF, p0 + HALF + 1
    if lo < 0 or hi > len(fa[k]):
        continue
    s = np.frombuffer(fa[k][lo:hi].encode(), dtype=np.uint8)
    s = CODE[s]
    # STRAND CONVENTION, MEASURED not assumed (diag_win.py, 8/8 sites):
    #   strand == 1  ->  reference base at pos-1 is G  ->  the site is on the MINUS strand
    #   strand == 0  ->  reference base at pos-1 is C  ->  plus strand, take as-is
    # My first attempt had this backwards and the built-in check caught it: 0.264 agreement
    # (chance is 0.25) with the focal base never C. That is bug 1's family -- a strand
    # convention correct in the script that wrote it, guessed wrong in the one that reads it.
    if int(strand[i]) == 1:
        s = COMP[s][::-1]
    out[i] = s
    if i % 200000 == 0 and i:
        print(f"  {i:,}/{N:,}")

# --- the verification that makes this trustworthy
c81 = HALF - mid81
sub = out[:, c81:c81 + win81.shape[1]]
idx = np.random.default_rng(0).choice(N, 200000, replace=False)
ok = (sub[idx] == win81[idx])
agree = float(ok.mean())
percol = ok.mean(axis=0)
print(f"\n1kb window vs existing 81bp window, 200,000 rows x {win81.shape[1]} bases: "
      f"{agree:.6f} agreement")
print(f"  worst column {percol.min():.6f} at offset {int(np.argmin(percol)) - mid81:+d}")
print(f"  focal base is C in {float((out[:, HALF] == 1).mean()):.6f} of rows")
if agree < 0.999:
    raise SystemExit("*** 1kb and 81bp windows DISAGREE -- coordinate or strand error, "
                     "not an architecture question. Nothing written. ***")

np.savez_compressed(OUT, win1k=out, chrom=chrom, pos=pos, strand=strand, y=d["y"])
print(f"wrote {OUT}  ({os.path.getsize(OUT)/1e6:.0f} MB)")
