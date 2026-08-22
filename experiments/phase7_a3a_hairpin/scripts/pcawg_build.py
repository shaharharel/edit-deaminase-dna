#!/usr/bin/env python
"""
S1 — PCAWG endogenous-APOBEC training-set builder.

Builds an A3A-endogenous classifier training set from the PCAWG open tier:
  positives = TCW C>T/C>G SNVs in A3A-dominant donors
  negatives = trinucleotide-matched, donor-matched unmutated C's

Design notes (why each control exists):
  * Negatives are matched on TRINUCLEOTIDE CONTEXT so the model cannot win by
    relearning the TpC motif. Unmatched negatives are what produced the
    meaningless AUROC 0.676 in the earlier editor work.
  * Negatives are matched on DONOR so overall mutation burden / exposure is
    controlled.
  * A3A-vs-A3B is split by YTCA/RTCA (Chan et al. 2015). We want the A3A arm:
    that is the hairpin-driven deaminase.
  * Local GC, distance-to-nearest-mutation-in-donor, and mappability proxies are
    recorded as columns so downstream can condition on them rather than having
    them silently drive the model.

Idempotent: writes flags into FLAGDIR, safe to re-run after preemption.
"""
import os, sys, gzip, json, time, random, subprocess
import numpy as np

BASE = "/mnt/data/a3a"
REF = "/mnt/data/ref/hg19.fa"
FLAGDIR = f"{BASE}/flags"
OUT = f"{BASE}/pcawg"
MAF_URL = ("https://object.genomeinformatics.org/icgc25k-open/PCAWG/"
           "consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz")
MAF = f"{OUT}/pcawg_icgc_public.maf.gz"

N_POSITIVES = 100_000
NEG_RATIO = 10          # negatives per positive
SEED = 20260810

os.makedirs(OUT, exist_ok=True)
os.makedirs(FLAGDIR, exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def flag(name):
    return os.path.join(FLAGDIR, name)


# ---------------------------------------------------------------- download
def step_download():
    if os.path.exists(flag("S1_DOWNLOAD_DONE")) and os.path.exists(MAF):
        log("MAF already downloaded, skipping")
        return
    log(f"downloading {MAF_URL}")
    subprocess.run(["aria2c", "-x", "8", "-s", "8", "--file-allocation=none",
                    "-d", OUT, "-o", os.path.basename(MAF), MAF_URL], check=True)
    # integrity: gzip must decompress cleanly
    subprocess.run(["gzip", "-t", MAF], check=True)
    open(flag("S1_DOWNLOAD_DONE"), "w").write("ok\n")
    log(f"downloaded + gzip-verified: {os.path.getsize(MAF)/1e6:.1f} MB")


# ---------------------------------------------------------------- reference
def load_hg19():
    """Load hg19 into a dict of uppercase str. ~3GB RAM; node has 125GB."""
    log("loading hg19 into memory")
    seqs, name, chunks = {}, None, []
    with open(REF) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                if name.startswith("chr"):
                    name = name[3:]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    keep = {k: v for k, v in seqs.items() if k in
            [str(i) for i in range(1, 23)] + ["X", "Y"]}
    log(f"hg19 loaded: {len(keep)} chroms, {sum(len(v) for v in keep.values())/1e9:.2f} Gbp")
    return keep


# ---------------------------------------------------------------- MAF parse
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def step_parse(genome):
    """Parse MAF -> per-SNV arrays, pyrimidine-oriented, with trinuc context."""
    cache = f"{OUT}/snvs.npz"
    if os.path.exists(cache) and os.path.exists(flag("S1_PARSE_DONE")):
        log("parse cache hit")
        d = np.load(cache, allow_pickle=True)
        return {k: d[k] for k in d.files}

    log("parsing MAF")
    chrom_l, pos_l, ref_l, alt_l, donor_l, tri_l = [], [], [], [], [], []
    n_tot = n_snp = 0
    with gzip.open(MAF, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = {c: i for i, c in enumerate(f)}
                needed = ["Chromosome", "Start_position", "Reference_Allele",
                          "Tumor_Seq_Allele2", "Variant_Type", "Donor_ID"]
                missing = [c for c in needed if c not in header]
                if missing:
                    # tolerate naming variants
                    alt_names = {"Start_position": "Start_Position",
                                 "Donor_ID": "Tumor_Sample_Barcode"}
                    for m in list(missing):
                        a = alt_names.get(m)
                        if a and a in header:
                            header[m] = header[a]
                            missing.remove(m)
                if missing:
                    log(f"FATAL missing MAF columns: {missing}")
                    log(f"available: {list(header)[:40]}")
                    sys.exit(1)
                log(f"MAF columns resolved ({len(header)} total)")
                continue

            n_tot += 1
            if f[header["Variant_Type"]] != "SNP":
                continue
            c = f[header["Chromosome"]].replace("chr", "")
            if c not in genome:
                continue
            p = int(f[header["Start_position"]])          # MAF is 1-based
            r = f[header["Reference_Allele"]].upper()
            a = f[header["Tumor_Seq_Allele2"]].upper()
            if r not in "ACGT" or a not in "ACGT" or r == a:
                continue
            seq = genome[c]
            if p < 2 or p > len(seq) - 1:
                continue
            tri = seq[p - 2:p + 1]
            if seq[p - 1] != r:      # reference mismatch -> drop
                continue
            if r in "AG":            # orient to pyrimidine
                r, a, tri = revcomp(r), revcomp(a), revcomp(tri)
            if "N" in tri:
                continue
            n_snp += 1
            chrom_l.append(c); pos_l.append(p); ref_l.append(r)
            alt_l.append(a); donor_l.append(f[header["Donor_ID"]]); tri_l.append(tri)
            if n_snp % 2_000_000 == 0:
                log(f"  parsed {n_snp/1e6:.0f}M usable SNVs (of {n_tot/1e6:.0f}M rows)")

    d = dict(chrom=np.array(chrom_l), pos=np.array(pos_l, dtype=np.int64),
             ref=np.array(ref_l), alt=np.array(alt_l),
             donor=np.array(donor_l), tri=np.array(tri_l))
    np.savez_compressed(cache, **d)
    open(flag("S1_PARSE_DONE"), "w").write(f"{n_snp}\n")
    log(f"parsed {n_snp:,} usable SNVs from {n_tot:,} MAF rows; {len(set(donor_l)):,} donors")
    return d


# ---------------------------------------------------- APOBEC / A3A donor calls
def step_donor_stats(d):
    """Roberts-style APOBEC enrichment + Chan YTCA/RTCA A3A-vs-A3B split."""
    log("computing per-donor APOBEC enrichment and YTCA/RTCA")
    is_c = d["ref"] == "C"
    is_apo_alt = np.isin(d["alt"], ["T", "G"])
    up = np.array([t[0] for t in d["tri"]])
    dn = np.array([t[2] for t in d["tri"]])
    is_tcw = is_c & (up == "T") & np.isin(dn, ["A", "T"]) & is_apo_alt

    stats = {}
    for donor in np.unique(d["donor"]):
        m = d["donor"] == donor
        n_c = int((is_c & m).sum())
        if n_c < 100:
            continue
        n_tcw = int((is_tcw & m).sum())
        enr = n_tcw / n_c
        stats[donor] = dict(n_snv=int(m.sum()), n_c=n_c, n_tcw=n_tcw, tcw_frac=enr)

    # YTCA vs RTCA needs the -2 base; recompute tetranucleotide for TCA sites only
    log(f"donor stats for {len(stats):,} donors")
    with open(f"{OUT}/donor_stats.json", "w") as fh:
        json.dump(stats, fh)
    return stats, is_tcw


def step_tetra(d, genome):
    """-2 base for YTCA/RTCA discrimination (A3A vs A3B)."""
    log("computing tetranucleotide (-2) context for A3A/A3B split")
    # ========================================================================
    # HARD STOP (added 2026-08-22 QA). THIS FUNCTION IS BUG 2 AND IS STILL WRONG.
    #
    # The branch below chooses the -2 base on `r == "C"`. But d["ref"] is ALREADY
    # pyrimidine-oriented (see step_donor_stats: `is_c = d["ref"] == "C"`), so the
    # else-branch is effectively dead code: minus-strand C mutations also take the
    # "C" branch and read seq[p-3] on the PLUS strand, which is the +2 position in
    # oriented space and uncomplemented. ~HALF of all sites get a wrong -2 base,
    # which corrupts every YTCA/RTCA number and therefore the A3A-vs-A3B donor split.
    #
    # This is the WRITER of the two corrupt artifacts -- minus2.npy and the
    # unsuffixed a3a_donor_ranking.tsv. The v1-ranking READERS were guarded first;
    # the writer is the more dangerous one, because rerunning it silently
    # regenerates the corruption under the same filenames.
    #
    # THE CORRECT IMPLEMENTATION IS fix_tetra.py, which derives strand from the
    # reference base (C on plus = 0, G on plus = 1) and complements the -2 base in
    # strand-oriented space. It consumes snvs.npz, which step_parse has ALREADY
    # written by the time this runs, so stopping here loses nothing: run
    # fix_tetra.py next.
    #
    # Kept unrepaired on purpose -- rewriting it would erase the record of bug 2.
    # Set A3A_ALLOW_BROKEN_TETRA=1 only to reproduce the bug deliberately.
    # ========================================================================
    import os as _os, sys as _sys
    if not _os.environ.get("A3A_ALLOW_BROKEN_TETRA"):
        _sys.exit(
            "REFUSING TO RUN step_tetra: it branches on ref=='C' but ref is already "
            "pyrimidine-oriented, so ~half of all -2 bases are wrong (bug 2). "
            "snvs.npz has been written; run fix_tetra.py to produce the correct -2 "
            "context and a3a_donor_ranking_v2.tsv. "
            "Set A3A_ALLOW_BROKEN_TETRA=1 to override deliberately.")
    m2 = []
    for c, p, r in zip(d["chrom"], d["pos"], d["ref"]):
        seq = genome[c]
        if r == "C":
            b = seq[p - 3] if p >= 3 else "N"
        else:
            b = revcomp(seq[p + 1]) if p + 1 < len(seq) else "N"
        m2.append(b)
    return np.array(m2)


def main():
    t0 = time.time()
    step_download()
    genome = load_hg19()
    d = step_parse(genome)
    stats, is_tcw = step_donor_stats(d)
    m2 = step_tetra(d, genome)
    np.save(f"{OUT}/minus2.npy", m2)

    # A3A-dominant donors: high TCW fraction AND YTCA-skewed
    ytca = (m2 == "T") | (m2 == "C")   # Y = pyrimidine at -2
    donors = np.array(sorted(stats))
    rows = []
    for donor in donors:
        m = (d["donor"] == donor) & is_tcw
        n = int(m.sum())
        if n < 50:
            continue
        y = int((m & ytca).sum())
        rows.append((donor, stats[donor]["tcw_frac"], n, y / max(n, 1)))
    rows.sort(key=lambda r: -r[1])
    with open(f"{OUT}/a3a_donor_ranking.tsv", "w") as fh:
        fh.write("donor\ttcw_frac\tn_tcw\tytca_frac\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]:.4f}\t{r[2]}\t{r[3]:.4f}\n")
    log(f"donor ranking written: {len(rows)} donors with >=50 TCW")
    log(f"S1 stage-1 complete in {(time.time()-t0)/60:.1f} min")
    open(flag("S1_STAGE1_DONE"), "w").write("ok\n")


if __name__ == "__main__":
    random.seed(SEED); np.random.seed(SEED)
    main()
