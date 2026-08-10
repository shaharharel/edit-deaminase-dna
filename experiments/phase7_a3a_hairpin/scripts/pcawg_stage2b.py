#!/usr/bin/env python
"""
S1 stage 2b — CORRECTED training-set build (strand-aware).

BUG FIXED: positives are pyrimidine-oriented (minus-strand mutations were
reverse-complemented), but stage2 sampled negatives by matching the trinucleotide
on the PLUS strand only. So every negative literally read TCA/TCT while ~half the
positives sat on a plus-strand genome reading TGA/AGA. Since revcomp(TCA)=TGA has
flanks T..A (same as TCA) but revcomp(TCT)=AGA has flanks A..A (vs T..T), a
sequence model could separate minus-strand TCT positives perfectly -- producing
exactly-ceiling 11.000x enrichment on a 100%-TCT top slice.

FIX: recover each positive's strand from the reference, sample negatives per
(context, strand) cell in the same proportion, and compute ALL features in
strand-oriented space (reverse-complement the window for minus-strand sites).
Positives and negatives are then exchangeable except for biology.
"""
import os, time
import numpy as np

BASE = "/mnt/data/a3a"
OUT = f"{BASE}/pcawg"
FEAT = f"{BASE}/feat"
REF = "/mnt/data/ref/hg19.fa"
N_POSITIVES, NEG_RATIO = 100_000, 10
TCW_MIN, YTCA_MIN = 0.20, 0.55
STEMS, LOOPS, FLANK = range(3, 11), range(3, 10), 40
SEED = 20260810
rng = np.random.default_rng(SEED)
B = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
COMPI = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def rcs(s):
    return "".join(RC[c] for c in reversed(s))


def load_hg19():
    log("loading hg19")
    seqs, name, chunks = {}, None, []
    with open(REF) as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0].replace("chr", "")
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks).upper()
    out = {}
    for c in [str(i) for i in range(1, 23)] + ["X"]:
        if c in seqs:
            a = np.frombuffer(seqs[c].encode(), dtype=np.uint8)
            o = np.full(a.shape, 4, dtype=np.uint8)
            for ch, v in B.items():
                if ch != "N":
                    o[a == ord(ch)] = v
            out[c] = o
    return out


def windows(genome, chrom, zpos, strand, flank):
    """Strand-oriented window of width 2*flank+1 centred on the focal base."""
    n = len(zpos)
    w = np.full((n, 2 * flank + 1), 4, dtype=np.uint8)
    for c in np.unique(chrom):
        m = chrom == c
        seq = genome[c]
        idx = np.clip(zpos[m][:, None] + np.arange(-flank, flank + 1)[None, :],
                      0, len(seq) - 1)
        ww = seq[idx]
        neg = strand[m] == 1
        if neg.any():
            ww[neg] = COMPI[ww[neg][:, ::-1]]
        w[m] = ww
    return w


def hairpin_from_windows(win):
    n = win.shape[0]
    C = win.shape[1] // 2
    best = np.full(n, -1.0)
    stem = np.zeros(n, np.int8); loop = np.zeros(n, np.int8)
    kpos = np.zeros(n, np.int8); gcp = np.zeros(n, np.int8)
    for L in LOOPS:
        for k in range(L):
            ls = C - k; le = ls + L
            for S in STEMS:
                a0, a1, b0, b1 = ls - S, ls, le, le + S
                if a0 < 0 or b1 > win.shape[1]:
                    continue
                left = win[:, a0:a1]; right = win[:, b0:b1]
                ok = (left == COMPI[right[:, ::-1]]).all(axis=1)
                if not ok.any():
                    continue
                gc = ((left == 1) | (left == 2)).sum(axis=1)
                sc = 2.0 * gc + 1.0 * (S - gc) - 0.25 * L
                u = ok & (sc > best)
                best[u] = sc[u]; stem[u] = S; loop[u] = L; kpos[u] = k; gcp[u] = gc[u]
    return stem, loop, kpos, gcp, best


def main():
    t0 = time.time()
    genome = load_hg19()
    d = np.load(f"{OUT}/snvs.npz", allow_pickle=True)
    chrom, pos, ref, alt = d["chrom"], d["pos"], d["ref"], d["alt"]
    donor, tri = d["donor"], d["tri"]

    rows = [l.split("\t") for l in
            open(f"{OUT}/a3a_donor_ranking.tsv").read().splitlines()[1:]]
    sel = {r[0] for r in rows if float(r[1]) >= TCW_MIN and float(r[3]) >= YTCA_MIN}
    in_sel = np.array([x in sel for x in donor])
    up = np.array([t[0] for t in tri]); dn = np.array([t[2] for t in tri])
    is_tcw = (ref == "C") & (up == "T") & np.isin(dn, ["A", "T"]) & np.isin(alt, ["T", "G"])
    cand = in_sel & is_tcw & np.isin(chrom, list(genome))
    idx = np.flatnonzero(cand)
    if len(idx) > N_POSITIVES:
        idx = rng.choice(idx, N_POSITIVES, replace=False)
    idx.sort()
    p_chrom, p_pos, p_donor = chrom[idx], pos[idx], donor[idx]
    p_tri = np.array([t[0] + t[1] + t[2] for t in tri[idx]])

    # ---- recover strand from the reference (C on plus = 0, G on plus = 1)
    p_strand = np.zeros(len(idx), np.int8)
    for c in np.unique(p_chrom):
        m = p_chrom == c
        base = genome[c][p_pos[m] - 1]
        p_strand[m] = np.where(base == 1, 0, 1)
    log(f"positives={len(idx):,}  plus-strand={int((p_strand==0).sum()):,} "
        f"minus-strand={int((p_strand==1).sum()):,}")
    bad = 0
    for c in np.unique(p_chrom):
        m = p_chrom == c
        base = genome[c][p_pos[m] - 1]
        bad += int(((base != 1) & (base != 2)).sum())
    log(f"  sanity: positives whose ref base is neither C nor G: {bad}")

    # ---- mutated positions per donor, to exclude from negatives
    mut = {}
    for c, p, dn_ in zip(chrom, pos, donor):
        mut.setdefault(dn_, set()).add((c, p))

    # ---- negatives matched on (trinucleotide context, STRAND) and donor
    log("sampling strand-aware trinuc-matched negatives")
    chrom_list = list(genome); clen = {c: len(genome[c]) for c in chrom_list}
    n_chrom, n_pos, n_donor, n_tri, n_strand = [], [], [], [], []
    for ctx in np.unique(p_tri):
        for st in (0, 1):
            cell = (p_tri == ctx) & (p_strand == st)
            need = int(cell.sum()) * NEG_RATIO
            if need == 0:
                continue
            donors_cell = p_donor[cell]
            want_s = ctx if st == 0 else rcs(ctx)
            want = np.array([B[want_s[0]], B[want_s[1]], B[want_s[2]]], np.uint8)
            got, tries = 0, 0
            while got < need and tries < 80:
                tries += 1
                batch = max(need * 25, 200_000)
                cc = rng.choice(chrom_list, batch)
                pp = np.zeros(batch, np.int64)
                for c in np.unique(cc):
                    m = cc == c
                    pp[m] = rng.integers(FLANK + 2, clen[c] - FLANK - 2, m.sum())
                keep = np.zeros(batch, bool)
                for c in np.unique(cc):
                    m = cc == c
                    s = genome[c]
                    keep[m] = ((s[pp[m] - 2] == want[0]) & (s[pp[m] - 1] == want[1])
                               & (s[pp[m]] == want[2]))
                for i in np.flatnonzero(keep):
                    if got >= need:
                        break
                    dn_ = donors_cell[got % len(donors_cell)]
                    if (cc[i], int(pp[i])) in mut.get(dn_, ()):
                        continue
                    n_chrom.append(cc[i]); n_pos.append(int(pp[i]))
                    n_donor.append(dn_); n_tri.append(ctx); n_strand.append(st)
                    got += 1
            log(f"  ctx={ctx} strand={st}: {got:,}/{need:,}")

    all_chrom = np.concatenate([p_chrom, np.array(n_chrom)])
    all_pos = np.concatenate([p_pos, np.array(n_pos)])
    all_donor = np.concatenate([p_donor, np.array(n_donor)])
    all_tri = np.concatenate([p_tri, np.array(n_tri)])
    all_str = np.concatenate([p_strand, np.array(n_strand, np.int8)])
    y = np.concatenate([np.ones(len(p_pos), np.int8), np.zeros(len(n_pos), np.int8)])
    zpos = all_pos - 1

    log("computing strand-oriented windows + hairpin features")
    win = windows(genome, all_chrom, zpos, all_str, FLANK)

    # LEAK CHECK: focal base must now be C for EVERY row, pos and neg alike
    focal = win[:, FLANK]
    log(f"  LEAK CHECK focal base==C: pos={float((focal[y==1]==1).mean()):.4f} "
        f"neg={float((focal[y==0]==1).mean()):.4f}  (both must be 1.0000)")
    for o in (-1, 1):
        fp = win[y == 1, FLANK + o]; fn = win[y == 0, FLANK + o]
        d1 = np.bincount(fp, minlength=5) / len(fp)
        d0 = np.bincount(fn, minlength=5) / len(fn)
        log(f"  LEAK CHECK offset {o:+d} ACGT pos={np.round(d1[:4],4)} "
            f"neg={np.round(d0[:4],4)}  (must match)")

    stem, loop, kpos, gcp, score = hairpin_from_windows(win)
    lgc = ((win == 1) | (win == 2)).mean(axis=1).astype(np.float32)

    np.savez_compressed(f"{FEAT}/a3a_trainset_v2.npz",
                        chrom=all_chrom, pos=all_pos, donor=all_donor, tri=all_tri,
                        strand=all_str, y=y, stem=stem, loop=loop, kpos=kpos,
                        gc_pairs=gcp, hp_score=score, local_gc=lgc, win=win)

    log("=" * 60)
    for s in range(3, 11):
        m = stem >= s
        p, n = m[y == 1].mean(), m[y == 0].mean()
        ys = rng.permutation(y)
        r = m[ys == 1].mean() / max(m[ys == 0].mean(), 1e-12)
        log(f"  stem>={s}: enr={p/max(n,1e-12):6.3f}x  RANDOM={r:5.3f}x  n_pos={int(m[y==1].sum()):,}")
    for k in range(6):
        m = (stem >= 4) & (kpos == k)
        p, n = m[y == 1].mean(), m[y == 0].mean()
        log(f"  stem>=4 & pos_in_loop={k}: enr={p/max(n,1e-12):6.3f}x")
    log("=" * 60)
    open(f"{BASE}/flags/S1_STAGE2B_DONE", "w").write("ok\n")
    log(f"stage2b complete in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
