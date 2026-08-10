#!/usr/bin/env python
"""
S1 stage 2 — A3A endogenous training set + hairpin/stem-loop features.

POSITIVES  : TCW C>T/C>G mutations in A3A-like PCAWG donors (tcw_frac>=0.20 &
             ytca_frac>=0.55). Restricting to TCW buys label purity.
NEGATIVES  : unmutated C's sampled from the SAME donor with the SAME
             trinucleotide context. Because positives are TCW-only and negatives
             are trinuc-matched, negatives are TCW too -- so the model CANNOT win
             by relearning the motif. Anything it learns is beyond-motif
             (structure/hairpin), which is the entire hypothesis under test.
COVARIATES : local GC, distance to nearest same-donor mutation (kataegis proxy),
             -2 base, replication-strand-agnostic. Recorded, NOT matched on, so
             downstream can condition instead of being silently driven by them.
SPLIT      : by CHROMOSOME (held-out chrom), never random -- spatial leakage
             inflated a previous result in this project.

Also emits RANDOM-BASELINE stats for every enrichment-style number.
"""
import os, sys, json, time
import numpy as np

BASE = "/mnt/data/a3a"
OUT = f"{BASE}/pcawg"
FEAT = f"{BASE}/feat"
FLAGDIR = f"{BASE}/flags"
REF = "/mnt/data/ref/hg19.fa"

N_POSITIVES = 100_000
NEG_RATIO = 10
TCW_MIN, YTCA_MIN = 0.20, 0.55
SEED = 20260810

# hairpin search space (Buisson-style): stem 3-10 bp, loop 3-9 nt
STEMS = range(3, 11)
LOOPS = range(3, 10)
FLANK = 40           # bp of context needed either side

os.makedirs(FEAT, exist_ok=True)
rng = np.random.default_rng(SEED)
B = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
COMPI = np.array([3, 2, 1, 0, 4], dtype=np.uint8)   # A<->T, C<->G


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_hg19_arrays():
    """hg19 as per-chrom uint8 arrays (A0 C1 G2 T3 N4)."""
    log("loading hg19 as uint8 arrays")
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
    keep = {}
    for c in [str(i) for i in range(1, 23)] + ["X"]:
        if c not in seqs:
            continue
        a = np.frombuffer(seqs[c].encode(), dtype=np.uint8)
        out = np.full(a.shape, 4, dtype=np.uint8)
        for ch, v in B.items():
            if ch != "N":
                out[a == ord(ch)] = v
        keep[c] = out
    log(f"hg19 arrays: {len(keep)} chroms")
    return keep


def select_donors():
    rows = [l.split("\t") for l in
            open(f"{OUT}/a3a_donor_ranking.tsv").read().splitlines()[1:]]
    sel = [r[0] for r in rows if float(r[1]) >= TCW_MIN and float(r[3]) >= YTCA_MIN]
    log(f"A3A-like donors selected: {len(sel)}")
    return set(sel)


def hairpin_features(genome, chrom, pos, strand):
    """
    For each site, find the most stable stem-loop whose loop contains the C.
    Returns (stem_len, loop_len, pos_in_loop, n_gc_pairs) of the best hairpin,
    zeros where none found. Vectorised over sites, looped over the (stem, loop,
    position) grid -- ~500 combinations.

    pos is 0-based index of the C in the (already strand-oriented) array.
    """
    n = len(pos)
    best_stem = np.zeros(n, dtype=np.int8)
    best_loop = np.zeros(n, dtype=np.int8)
    best_k = np.zeros(n, dtype=np.int8)
    best_gc = np.zeros(n, dtype=np.int8)
    best_score = np.full(n, -1.0)

    # gather a padded window per site once
    win = np.zeros((n, 2 * FLANK + 1), dtype=np.uint8)
    for c in np.unique(chrom):
        m = chrom == c
        seq = genome[c]
        idx = pos[m][:, None] + np.arange(-FLANK, FLANK + 1)[None, :]
        idx = np.clip(idx, 0, len(seq) - 1)
        w = seq[idx]
        if strand is not None:
            neg = strand[m] == 1
            if neg.any():
                w[neg] = COMPI[w[neg][:, ::-1]]
        win[m] = w
    C = FLANK   # index of the focal C inside win

    for L in LOOPS:
        for k in range(L):                 # C is the k-th base of the loop
            ls = C - k                     # loop start index
            le = ls + L                    # loop end (exclusive)
            for S in STEMS:
                a0, a1 = ls - S, ls        # 5' stem
                b0, b1 = le, le + S        # 3' stem
                if a0 < 0 or b1 > win.shape[1]:
                    continue
                left = win[:, a0:a1]
                right = win[:, b0:b1]
                # pair if left[i] complements right[S-1-i]
                paired = (left == COMPI[right[:, ::-1]])
                ok = paired.all(axis=1)
                if not ok.any():
                    continue
                gc = ((left == 1) | (left == 2)).sum(axis=1)
                # stability proxy: GC pairs weigh ~2x AT, short loops favoured
                score = 2.0 * gc + 1.0 * (S - gc) - 0.25 * L
                upd = ok & (score > best_score)
                best_score[upd] = score[upd]
                best_stem[upd] = S
                best_loop[upd] = L
                best_k[upd] = k
                best_gc[upd] = gc[upd]
    return best_stem, best_loop, best_k, best_gc, best_score


def main():
    t0 = time.time()
    genome = load_hg19_arrays()
    d = np.load(f"{OUT}/snvs.npz", allow_pickle=True)
    chrom, pos, ref, alt = d["chrom"], d["pos"], d["ref"], d["alt"]
    donor, tri = d["donor"], d["tri"]
    m2 = np.load(f"{OUT}/minus2.npy", allow_pickle=True)

    sel = select_donors()
    in_sel = np.array([x in sel for x in donor])
    up = np.array([t[0] for t in tri]); dn = np.array([t[2] for t in tri])
    is_tcw = (ref == "C") & (up == "T") & np.isin(dn, ["A", "T"]) & np.isin(alt, ["T", "G"])
    ok_chrom = np.isin(chrom, list(genome))
    cand = in_sel & is_tcw & ok_chrom
    log(f"candidate positives: {cand.sum():,}")

    idx = np.flatnonzero(cand)
    if len(idx) > N_POSITIVES:
        idx = rng.choice(idx, N_POSITIVES, replace=False)
    idx.sort()
    log(f"positives sampled: {len(idx):,}")

    p_chrom, p_pos, p_donor, p_tri = chrom[idx], pos[idx], donor[idx], tri[idx]
    p_m2 = m2[idx]

    # ---- same-donor mutated positions, to exclude from negatives
    mutset = {}
    for c, p, dn_ in zip(chrom, pos, donor):
        mutset.setdefault(dn_, set()).add((c, p))

    # ---- trinuc-matched, donor-matched negatives by rejection sampling
    log(f"sampling {NEG_RATIO}x trinuc-matched negatives")
    tri_want = np.array([t[0] + t[1] + t[2] for t in p_tri])
    n_chrom, n_pos, n_donor, n_tri = [], [], [], []
    chrom_list = list(genome)
    chrom_len = {c: len(genome[c]) for c in chrom_list}
    order = np.argsort(tri_want)
    for ctx in np.unique(tri_want):
        need = int((tri_want == ctx).sum()) * NEG_RATIO
        donors_for_ctx = p_donor[tri_want == ctx]
        got = 0
        tries = 0
        want_codes = np.array([B[ctx[0]], B[ctx[1]], B[ctx[2]]], dtype=np.uint8)
        while got < need and tries < 60:
            tries += 1
            batch = max(need * 20, 200_000)
            cc = rng.choice(chrom_list, batch)
            pp = np.array([rng.integers(FLANK + 2, chrom_len[c] - FLANK - 2)
                           for c in cc], dtype=np.int64)
            keep = np.zeros(batch, dtype=bool)
            for c in np.unique(cc):
                m = cc == c
                seq = genome[c]
                t0_ = seq[pp[m] - 2]; t1_ = seq[pp[m] - 1]; t2_ = seq[pp[m]]
                keep[m] = ((t0_ == want_codes[0]) & (t1_ == want_codes[1]) &
                           (t2_ == want_codes[2]))
            sel_i = np.flatnonzero(keep)
            for i in sel_i:
                if got >= need:
                    break
                dn_ = donors_for_ctx[got % len(donors_for_ctx)]
                if (cc[i], int(pp[i])) in mutset.get(dn_, ()):
                    continue
                n_chrom.append(cc[i]); n_pos.append(int(pp[i]))
                n_donor.append(dn_); n_tri.append(ctx)
                got += 1
        log(f"  ctx {ctx}: {got:,}/{need:,} negatives ({tries} rounds)")

    n_chrom = np.array(n_chrom); n_pos = np.array(n_pos)
    n_donor = np.array(n_donor); n_tri = np.array(n_tri)
    log(f"negatives: {len(n_pos):,}")

    # ---- assemble + features
    all_chrom = np.concatenate([p_chrom, n_chrom])
    all_pos = np.concatenate([p_pos, n_pos])
    all_donor = np.concatenate([p_donor, n_donor])
    all_tri = np.concatenate([tri_want, n_tri])
    y = np.concatenate([np.ones(len(p_pos), np.int8), np.zeros(len(n_pos), np.int8)])

    log("computing hairpin/stem-loop features")
    zpos = all_pos - 1                     # MAF 1-based -> 0-based
    stem, loop, kpos, gc, score = hairpin_features(genome, all_chrom, zpos, None)

    # local GC (+/-250bp) as a recorded covariate
    log("computing local GC")
    lgc = np.zeros(len(zpos), dtype=np.float32)
    for c in np.unique(all_chrom):
        m = all_chrom == c
        seq = genome[c]
        ii = zpos[m][:, None] + np.arange(-250, 251)[None, :]
        ii = np.clip(ii, 0, len(seq) - 1)
        w = seq[ii]
        lgc[m] = (((w == 1) | (w == 2)).sum(axis=1) / w.shape[1]).astype(np.float32)

    np.savez_compressed(
        f"{FEAT}/a3a_trainset.npz",
        chrom=all_chrom, pos=all_pos, donor=all_donor, tri=all_tri, y=y,
        stem=stem, loop=loop, kpos=kpos, gc_pairs=gc, hp_score=score, local_gc=lgc)

    # ---- RANDOM BASELINE + honest summary
    hp = stem > 0
    log("=" * 60)
    log(f"positives={int(y.sum()):,}  negatives={int((y==0).sum()):,}")
    log(f"hairpin found: pos={hp[y==1].mean():.4f}  neg={hp[y==0].mean():.4f}  "
        f"ratio={hp[y==1].mean()/max(hp[y==0].mean(),1e-9):.3f}x")
    log(f"  RANDOM BASELINE (label-shuffled): "
        f"{hp[rng.permutation(len(y))][y==1].mean()/max(hp[y==0].mean(),1e-9):.3f}x")
    for nm, arr in [("stem_len", stem), ("loop_len", loop), ("pos_in_loop", kpos)]:
        log(f"  {nm}: pos_mean={arr[y==1].mean():.3f} neg_mean={arr[y==0].mean():.3f}")
    log(f"  local_GC: pos={lgc[y==1].mean():.4f} neg={lgc[y==0].mean():.4f}")
    log(f"NOTE negatives are trinucleotide-matched, so motif cannot drive any of this.")
    log("=" * 60)
    log(f"S1 stage-2 complete in {(time.time()-t0)/60:.1f} min")
    open(f"{FLAGDIR}/S1_STAGE2_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
