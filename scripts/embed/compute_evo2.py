"""Compute Evo-2 embeddings for DNA off-target sites on A100.

Evo-2 1B is the DNA-RNA multi-modal foundation model from Arc Institute (Brixi 2025).
Requires compute capability ≥ 8.0 (A100 / H100).

Inputs (on A100-b at ~/dna_emb/):
  - dna_offtarget_sites.parquet (53K sites)
  - hg38.fa (3.1 GB, pyfaidx-indexed)

Output:
  - evo2_dna_offtarget.pt  dict[site_key -> torch.Tensor]

Embedding strategy: fetch ±512 bp around (chrom, pos), tokenize via Evo-2's
byte/char tokenizer, forward through model, mean-pool last hidden states.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import pysam

WINDOW = 512
MODEL_ID = "arcinstitute/evo2_1b_base"


def get_context(fa: pysam.FastaFile, chrom: str, pos: int, window: int = WINDOW) -> str:
    try:
        start = max(0, pos - 1 - window)
        end = pos - 1 + window + 1
        return fa.fetch(chrom, start, end).upper()
    except (KeyError, ValueError):
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="~/dna_emb/dna_offtarget_sites.parquet")
    ap.add_argument("--hg38", default="~/dna_emb/hg38.fa")
    ap.add_argument("--out", default="~/dna_emb/evo2_dna_offtarget.pt")
    ap.add_argument("--batch", type=int, default=4,
                    help="Evo-2 is large — small batches on 40 GB A100")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    pq = Path(args.parquet).expanduser()
    fa_path = Path(args.hg38).expanduser()
    out = Path(args.out).expanduser()

    df = pd.read_parquet(pq)
    if args.limit > 0:
        df = df.head(args.limit)
    valid = df["chrom"].str.match(r"^chr(\d+|[XYM])$")
    df = df[valid].reset_index(drop=True)
    print(f"[evo2] {len(df)} sites")

    fa = pysam.FastaFile(str(fa_path))
    seqs, keys = [], []
    for _, r in df.iterrows():
        seq = get_context(fa, r["chrom"], int(r["pos"]), WINDOW)
        if len(seq) < 200:
            continue
        keys.append(f"{r['chrom']}:{int(r['pos'])}:{r['strand']}:{r['ref']}:{r['alt']}")
        seqs.append(seq)
    print(f"[evo2] {len(seqs)} sequences fetched (mean len {np.mean([len(s) for s in seqs]):.0f})")

    # Try the official Evo-2 wrapper first (cleanest API)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[evo2] device={device}, capability={torch.cuda.get_device_capability(0)}")

    try:
        from evo2 import Evo2
        print(f"[evo2] using Evo2 official wrapper")
        model = Evo2("evo2_1b_base")
        # Evo2 wrapper exposes .embed() or .forward(seq, return_hidden=True)
        # Embed in singleton batches first to be safe on memory
        embeddings = {}
        for i, (k, s) in enumerate(zip(keys, seqs)):
            with torch.no_grad():
                # Evo2 wrapper API: .embed(seq) returns per-position embedding [L, D]
                emb_per_pos = model.embed(s)
                emb = emb_per_pos.mean(dim=0).cpu().float()
            embeddings[k] = emb
            if i % 50 == 0:
                print(f"  {i+1}/{len(seqs)}")
    except ImportError:
        print(f"[evo2] evo2 wrapper unavailable — falling back to transformers AutoModel")
        from transformers import AutoModelForCausalLM, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
        model = AutoModelForCausalLM.from_pretrained(
            MODEL_ID, trust_remote_code=True, torch_dtype=torch.bfloat16).to(device).eval()
        n = len(seqs)
        embeddings = {}
        for i in range(0, n, args.batch):
            bk = keys[i:i + args.batch]
            bs = seqs[i:i + args.batch]
            enc = tok(bs, return_tensors="pt", padding=True, truncation=True, max_length=1100)
            enc = {k: v.to(device) for k, v in enc.items()}
            with torch.no_grad():
                out_h = model(**enc, output_hidden_states=True)
            h = out_h.hidden_states[-1]
            emb = h.mean(dim=1).detach().cpu().to(torch.float32)
            for k, e in zip(bk, emb):
                embeddings[k] = e.clone()
            if (i // args.batch) % 20 == 0:
                print(f"  {i+len(bk)}/{n}")

    torch.save(embeddings, out)
    sample_dim = next(iter(embeddings.values())).shape
    print(f"[done] wrote {out} — {len(embeddings)} embeddings, dim {sample_dim}")


if __name__ == "__main__":
    main()
