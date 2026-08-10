#!/usr/bin/env python
"""
S6 — per-clone pileup over the fixed scoring universe.

Turns an aligned clone BAM into (coverage, alt_count) aligned index-for-index to
feat/universe_chr1.npz, so every clone is scored on the SAME reference-derived
site list. The universe was built before any clone existed, which is what keeps
the later enrichment honest.

COORDINATE SAFETY (see CONVENTIONS.md): universe `pos` is 0-based, `pos1` is
1-based. samtools emits 1-based. We join on pos1 and then ASSERT the reference
base is C (strand 0) or G (strand 1) at a sample of joined sites. An off-by-one
here does not crash -- it silently turns TCW into CWN and yields a plausible
wrong answer. That assertion is the whole point.

Strand-aware counting: plus-strand sites count C->T, minus-strand sites count
G->A. Counting only one of these is how a 21.8:1 strand skew got into an earlier
stream in this project.

Usage: s6_pileup.py <sample> <bam> [chrom]
"""
import os, sys, subprocess, time
import numpy as np
sys.path.insert(0, '/mnt/data/a3a')
from pileup_parse import count_alt, count_depth

BASE = "/mnt/data/a3a"
FEAT = f"{BASE}/feat"
REF = "/mnt/data/ref/hg19.fa"
SAMTOOLS = os.path.expanduser("~/miniconda3/envs/bio/bin/samtools")
MIN_MQ, MIN_BQ, MAX_DEPTH = 20, 20, 500

sample, bam = sys.argv[1], sys.argv[2]
chrom = sys.argv[3] if len(sys.argv) > 3 else "1"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {sample}: {m}", flush=True)



def assert_same_build(bam, ref):
    """Abort unless the BAM's @SQ lengths match the reference .fai.

    THIS is the check that catches a BAM aligned to the wrong genome. The
    per-site reference-base comparison does NOT: samtools reports the reference
    base from the -f FASTA, not from the reads, so comparing it to a universe
    built from the same FASTA is a tautology (verified: hg38 BAM vs hg19 universe
    gave 0/158,973 mismatches).
    """
    fai = {}
    for line in open(ref + ".fai"):
        f = line.split("\t")
        fai[f[0]] = int(f[1])
    hdr = subprocess.run([SAMTOOLS, "view", "-H", bam], capture_output=True,
                         text=True).stdout
    checked = bad = 0
    for line in hdr.splitlines():
        if not line.startswith("@SQ"):
            continue
        d = dict(kv.split(":", 1) for kv in line.split("\t")[1:] if ":" in kv)
        sn, ln = d.get("SN"), int(d.get("LN", -1))
        if sn in fai:
            checked += 1
            if fai[sn] != ln:
                bad += 1
                if bad <= 3:
                    log(f"  BUILD MISMATCH {sn}: bam={ln:,} ref={fai[sn]:,}")
    if checked == 0:
        log("FATAL: no @SQ names shared between BAM and reference")
        sys.exit(3)
    if bad:
        log(f"FATAL: {bad}/{checked} sequences differ in length -- BAM is aligned "
            f"to a DIFFERENT BUILD than {ref}. Refusing to pile up.")
        sys.exit(3)
    log(f"  build check OK: {checked} sequences match {ref}")

def main():
    out = f"{FEAT}/counts_{sample}_chr{chrom}.npz"
    if os.path.exists(out):
        log("counts exist, skip"); return

    assert_same_build(bam, REF)

    u = np.load(f"{FEAT}/universe_chr{chrom}.npz")
    pos1 = u["pos1"].astype(np.int64)      # 1-based, for the samtools join
    strand = u["strand"]
    log(f"universe: {len(pos1):,} sites")

    # index from 1-based coordinate -> universe row
    # moving pointer instead of a dict: universe is coordinate-sorted and
    # mpileup emits in coordinate order, so the join is O(n) with no 2 GB dict.
    ptr = 0
    n_pos = len(pos1)
    n_unmatched = 0

    cov = np.zeros(len(pos1), np.int32)
    alt = np.zeros(len(pos1), np.int32)
    alt_f = np.zeros(len(pos1), np.int32)
    alt_r = np.zeros(len(pos1), np.int32)

    cmd = [SAMTOOLS, "mpileup", "-f", REF, "-l", f"{FEAT}/universe_chr{chrom}.bed",
           "-r", f"chr{chrom}", "-q", str(MIN_MQ), "-Q", str(MIN_BQ),
           "-d", str(MAX_DEPTH), "--no-BAQ", bam]
    log("running mpileup")
    n_seen = n_join = 0
    checked = 0
    bad_ref = 0
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)
    for line in proc.stdout:
        f = line.split("\t")
        if len(f) < 5:
            continue
        p = int(f[1])                      # samtools: 1-based
        while ptr < n_pos and pos1[ptr] < p:
            ptr += 1
        if ptr >= n_pos:
            break
        if pos1[ptr] != p:
            n_unmatched += 1               # emitted a site not in the universe
            continue
        i = ptr
        n_seen += 1
        rb = f[2].upper()
        st = strand[i]
        # ---- COORDINATE ASSERTION (sampled, cheap)
        if checked < 5000:
            checked += 1
            want = "C" if st == 0 else "G"
            if rb != want:
                bad_ref += 1
        bases = f[4]
        # BOTH strands must be counted. Uppercase = forward-strand read,
        # lowercase = reverse-strand read; counting only one loses ~half the
        # support ASYMMETRICALLY -- the mechanism behind the 21.8:1 skew that
        # contaminated the earlier Selict stream. Indel and read-start tokens are
        # stripped by pileup_parse; raw .count() would score them as bases.
        want = "T" if st == 0 else "A"     # plus: C->T ; minus: G->A on plus
        fw, rv = count_alt(bases, want)
        cov[i] = count_depth(bases)
        alt[i] = fw + rv
        alt_f[i] = fw
        alt_r[i] = rv
        n_join += 1
    proc.wait()

    if n_unmatched > 0.01 * max(n_join, 1):
        log(f"WARNING: {n_unmatched:,} mpileup sites absent from the universe "
            f"({n_unmatched/max(n_join,1):.3f} of joined) -- BED/universe mismatch?")
    frac_bad = bad_ref / max(checked, 1)
    log(f"joined {n_join:,} sites; universe-vs-reference agreement "
        f"(NOT a BAM check -- see assert_same_build): {bad_ref}/{checked} bad "
        f"({frac_bad:.4f})")
    if frac_bad > 0.01:
        log("FATAL: reference base mismatch > 1% -- coordinate convention is wrong. "
            "Refusing to write counts (see CONVENTIONS.md).")
        sys.exit(2)

    covered = cov >= 8
    log(f"sites cov>=8: {int(covered.sum()):,} ({covered.mean():.3f})")
    log(f"sites with alt>=1: {int((alt>=1).sum()):,}; alt>=2: {int((alt>=2).sum()):,}")
    tf, tr = int(alt_f.sum()), int(alt_r.sum())
    ratio = tf / max(tr, 1)
    log(f"STRAND BALANCE alt fwd={tf:,} rev={tr:,} ratio={ratio:.3f} "
        f"(expect ~1.0; a large skew means a counting/orientation bug)")
    if tf + tr > 1000 and (ratio > 2.0 or ratio < 0.5):
        log("WARNING: strand-asymmetric alt counts -- investigate before use.")
    np.savez_compressed(out, cov=cov, alt=alt, alt_fwd=alt_f, alt_rev=alt_r,
                        sample=np.array([sample]))
    open(f"{BASE}/flags/S6_{sample}_chr{chrom}_DONE", "w").write("ok\n")
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
