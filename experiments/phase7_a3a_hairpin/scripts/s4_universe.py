#!/usr/bin/env python
"""
S4a — build the scoring universe: every TCW cytosine (both strands) on the
target chromosomes, with strand-oriented hairpin features precomputed.

This is the denominator for the editor test. Defining it ONCE, independent of any
clone's calls, is what keeps the later enrichment honest: sites are ranked by
model score computed from the reference alone, never from the data being scored.

Output: feat/universe_<chr>.npz  (pos, strand, tri, stem, loop, kpos, gc_pairs,
hp_score, local_gc, win) + a BED for samtools mpileup.
"""
import os, sys, time
import numpy as np

BASE = "/mnt/data/a3a"
FEAT = f"{BASE}/feat"
REF = "/mnt/data/ref/hg19.fa"
CHROMS = sys.argv[1:] or ["1"]
STEMS, LOOPS, FLANK = range(3, 11), range(3, 10), 40
COMPI = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
B = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_chrom(c):
    want = f"chr{c}"
    seq, name, chunks = None, None, []
    with open(REF) as fh:
        for line in fh:
            if line[0] == ">":
                if name in (want, c) and chunks:
                    seq = "".join(chunks).upper(); break
                name = line[1:].split()[0]
                chunks = []
            elif name in (want, c):
                chunks.append(line.strip())
    if seq is None and chunks:
        seq = "".join(chunks).upper()
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    o = np.full(a.shape, 4, dtype=np.uint8)
    for ch, v in B.items():
        if ch != "N":
            o[a == ord(ch)] = v
    return o


def hairpin(win):
    n, W = win.shape
    C = W // 2
    best = np.full(n, -1.0, dtype=np.float32)
    stem = np.zeros(n, np.int8); loop = np.zeros(n, np.int8)
    kpos = np.zeros(n, np.int8); gcp = np.zeros(n, np.int8)
    for L in LOOPS:
        for k in range(L):
            ls = C - k; le = ls + L
            for S in STEMS:
                a0, a1, b0, b1 = ls - S, ls, le, le + S
                if a0 < 0 or b1 > W:
                    continue
                left = win[:, a0:a1]; right = win[:, b0:b1]
                ok = (left == COMPI[right[:, ::-1]]).all(axis=1)
                if not ok.any():
                    continue
                gc = ((left == 1) | (left == 2)).sum(axis=1)
                sc = (2.0 * gc + 1.0 * (S - gc) - 0.25 * L).astype(np.float32)
                u = ok & (sc > best)
                best[u] = sc[u]; stem[u] = S; loop[u] = L; kpos[u] = k; gcp[u] = gc[u]
    return stem, loop, kpos, gcp, best


def main():
    for c in CHROMS:
        out = f"{FEAT}/universe_chr{c}.npz"
        if os.path.exists(out):
            log(f"chr{c} universe exists, skip"); continue
        t0 = time.time()
        seq = load_chrom(c)
        n = len(seq)
        log(f"chr{c}: {n/1e6:.1f} Mb")

        # TCW on plus strand: T C [A|T]  -> focal C at i
        i = np.arange(FLANK + 1, n - FLANK - 1)
        plus = (seq[i - 1] == 3) & (seq[i] == 1) & ((seq[i + 1] == 0) | (seq[i + 1] == 3))
        # TCW on minus strand: plus reads [A|T] G A  -> focal G at i
        minus = (seq[i + 1] == 0) & (seq[i] == 2) & ((seq[i - 1] == 3) | (seq[i - 1] == 0))
        pos = np.concatenate([i[plus], i[minus]])
        strand = np.concatenate([np.zeros(plus.sum(), np.int8),
                                 np.ones(minus.sum(), np.int8)])
        o = np.argsort(pos); pos, strand = pos[o], strand[o]
        log(f"  TCW sites: {len(pos):,} (plus {int((strand==0).sum()):,} / "
            f"minus {int((strand==1).sum()):,})")

        # strand-oriented windows, chunked to bound memory
        stem = np.zeros(len(pos), np.int8); loop = np.zeros(len(pos), np.int8)
        kpos = np.zeros(len(pos), np.int8); gcp = np.zeros(len(pos), np.int8)
        hps = np.zeros(len(pos), np.float32); lgc = np.zeros(len(pos), np.float32)
        tri = np.zeros(len(pos), np.int8)   # 0=TCA 1=TCT
        CH = 2_000_000
        for s0 in range(0, len(pos), CH):
            s1 = min(s0 + CH, len(pos))
            idx = pos[s0:s1][:, None] + np.arange(-FLANK, FLANK + 1)[None, :]
            w = seq[np.clip(idx, 0, n - 1)]
            neg = strand[s0:s1] == 1
            if neg.any():
                w[neg] = COMPI[w[neg][:, ::-1]]
            assert (w[:, FLANK] == 1).all(), "focal base must be C after orientation"
            tri[s0:s1] = (w[:, FLANK + 1] == 3).astype(np.int8)   # 0=A(TCA) 1=T(TCT)
            st, lo, kp, gp, sc = hairpin(w)
            stem[s0:s1], loop[s0:s1], kpos[s0:s1] = st, lo, kp
            gcp[s0:s1], hps[s0:s1] = gp, sc
            lgc[s0:s1] = ((w == 1) | (w == 2)).mean(axis=1)
            log(f"  featurised {s1:,}/{len(pos):,}")

        # pos is the 0-based array index; pos1 is 1-based for joins against
        # samtools/MAF. BOTH are written so downstream never has to guess -- a
        # fallback to `pos` at a 1-based join shifts every site by one base and
        # silently turns TCW into CWN (see CONVENTIONS.md).
        np.savez_compressed(out, pos=pos.astype(np.int32),
                            pos1=(pos.astype(np.int64) + 1), strand=strand, tri=tri,
                            stem=stem, loop=loop, kpos=kpos, gc_pairs=gcp,
                            hp_score=hps, local_gc=lgc,
                            coord_note=np.array(['pos=0-based; pos1=1-based']))
        with open(f"{FEAT}/universe_chr{c}.bed", "w") as fh:
            for p in pos:
                fh.write(f"chr{c}\t{p}\t{p+1}\n")
        log(f"chr{c} done in {(time.time()-t0)/60:.1f} min -> {out}")
        log(f"  hairpin stem>=6 fraction: {(stem>=6).mean():.5f}  "
            f"stem>=8: {(stem>=8).mean():.6f}")
    open(f"{BASE}/flags/S4_UNIVERSE_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
