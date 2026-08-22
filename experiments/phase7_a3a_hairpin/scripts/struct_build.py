#!/usr/bin/env python
"""Thermodynamic structure modality for DeaminaFormer-DNA.

Replaces a scanner's stem-length integer with the quantity the mechanism actually names:
A3A deaminates a cytosine only when that C is SINGLE STRANDED. The partition function gives
p_unpaired(focal) directly.

DNA PARAMETERS, measured not assumed. On the same 150 windows, Mathews2004 DNA gives mean
MFE -8.9 kcal/mol where Turner2004 RNA gives -21.1 -- a 2.4x difference. Folding ssDNA with
RNA energies would have been a quiet, plausible-looking error of exactly the family this
project keeps finding.

TWO SCALES, because the answer depends on the window: p_unpaired(focal) is 0.563 at +-50 bp
and 0.488 at +-100 bp. A single scale would hide that dependence; two make it a feature.
+-25 and +-50 are used -- APOBEC hairpin stems are <=11 bp and loops <=11 nt, so local
windows describe the substrate, while wide windows let distant structure dominate the MFE.
+-100 was measured at 35 core-hours and is not worth 3.4x the cost for a diluted signal.

PRE-REGISTERED, and weak: on a 600-site probe, positives were more unpaired than negatives
(0.5195 vs 0.4903, +0.097 pooled SD) -- the direction the ssDNA mechanism predicts, but on
only 65 positives. That is the n where this project's A3A-vs-A3B claim died. This block is
built because it is mechanistically right, NOT because a 65-positive probe was encouraging.
"""
import os, sys, time
import numpy as np
import RNA
from multiprocessing import Pool

BASE = "/mnt/a3a"
HALVES = (25, 50)
NPROC = int(os.environ.get("STRUCT_PROC", "6"))
LUT = np.array(list("ACGTN"))

W = np.load(f"{BASE}/a3a_trainset_v5_win1k.npz", allow_pickle=True)
win = W["win1k"]; N = len(W["y"]); MID = win.shape[1] // 2


def init():
    RNA.params_load_DNA_Mathews2004()


def one(i):
    row = win[i]
    out = []
    for H in HALVES:
        s = "".join(LUT[row[MID - H:MID + H + 1]])
        fc = RNA.fold_compound(s)
        st, e_mfe = fc.mfe()
        fc.exp_params_rescale(e_mfe)
        e_ens = fc.pf()[1]
        bpp = np.array(fc.bpp())
        j = H + 1
        pu = 1.0 - float(bpp[j, :].sum() + bpp[:, j].sum())
        lo, hi = max(j - 5, 1), min(j + 5, 2 * H + 1)
        pu5 = float(np.mean([1.0 - (bpp[k, :].sum() + bpp[:, k].sum()) for k in range(lo, hi + 1)]))
        out += [pu, pu5, e_mfe, e_ens, float(st[j - 1] == ".")]
    return out


if __name__ == "__main__":
    t0 = time.time()
    with Pool(NPROC, initializer=init) as p:
        res = []
        for k, r in enumerate(p.imap(one, range(N), chunksize=512)):
            res.append(r)
            if k % 50000 == 0 and k:
                el = time.time() - t0
                print(f"  {k:,}/{N:,}  {k/el:.0f} sites/s  ETA {(N-k)/(k/el)/60:.0f} min",
                      flush=True)
    A = np.asarray(res, dtype=np.float32)
    names = [f"{n}_h{H}" for H in HALVES for n in ("p_unpaired", "p_unpaired5", "mfe", "e_ens", "mfe_unpaired")]
    print(f"structure block {A.shape} in {(time.time()-t0)/60:.1f} min")
    for c, n in enumerate(names):
        print(f"  {n:22s} mean {A[:, c].mean():8.3f}  sd {A[:, c].std():7.3f}")
    np.savez(f"{BASE}/struct_thermo_v5.npz", X=A, names=np.array(names),
             chrom=W["chrom"], pos=W["pos"], y=W["y"])
    print(f"wrote {BASE}/struct_thermo_v5.npz")
