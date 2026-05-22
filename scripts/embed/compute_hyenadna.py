"""Compute HyenaDNA embeddings for DNA off-target sites (DeaminaFormer-DNA features).

Runs on T4 (~/dna_emb/). Reads dna_offtarget_sites.parquet, extracts a ±256 bp window
from hg38 around each (chrom, pos), embeds with HyenaDNA-small (3.3M params, 256d).

Outputs: hyenadna_dna_offtarget.pt
    dict[site_key (chrom:pos:strand:ref:alt)] -> torch.Tensor[256]

Plan B if Evo unavailable on T4. Smaller than Evo embeddings but task-suited.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import pysam
from transformers import AutoModel, AutoTokenizer

MODEL_ID = "LongSafari/hyenadna-small-32k-seqlen-hf"
WINDOW = 256  # ±256 bp = 513 bp context


def get_context(fa: pysam.FastaFile, chrom: str, pos: int, window: int = WINDOW) -> str:
    try:
        # pysam is 0-based; pos in parquet appears 1-based (SRA standard)
        start = max(0, pos - 1 - window)
        end = pos - 1 + window + 1
        return fa.fetch(chrom, start, end).upper()
    except (KeyError, ValueError):
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="~/dna_emb/dna_offtarget_sites.parquet")
    ap.add_argument("--hg38", default="~/dna_emb/hg38.fa")
    ap.add_argument("--out", default="~/dna_emb/hyenadna_dna_offtarget.pt")
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--limit", type=int, default=0, help="0 = all sites")
    args = ap.parse_args()

    pq = Path(args.parquet).expanduser()
    fa_path = Path(args.hg38).expanduser()
    out = Path(args.out).expanduser()

    df = pd.read_parquet(pq)
    if args.limit > 0:
        df = df.head(args.limit)
    print(f"[embed] {len(df)} sites total")

    # Drop non-standard chroms
    valid = df["chrom"].str.match(r"^chr(\d+|[XYM])$")
    df = df[valid].reset_index(drop=True)
    print(f"[embed] {len(df)} sites after chrom filter")

    # Build the lookup key + fetch sequences (single-thread; pysam is fast)
    fa = pysam.FastaFile(str(fa_path))
    print(f"[embed] FASTA refs: {len(fa.references)}")
    seqs = []
    keys = []
    for _, r in df.iterrows():
        chrom = r["chrom"]; pos = int(r["pos"])
        seq = get_context(fa, chrom, pos, WINDOW)
        if len(seq) < 100:
            continue
        keys.append(f"{chrom}:{pos}:{r['strand']}:{r['ref']}:{r['alt']}")
        seqs.append(seq)
    print(f"[embed] fetched {len(seqs)} sequences (mean len {np.mean([len(s) for s in seqs]):.0f})")

    # Load model
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[embed] device={device}, loading {MODEL_ID}")
    tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
    model = AutoModel.from_pretrained(MODEL_ID, trust_remote_code=True).to(device).eval()
    print(f"[embed] params (M): {sum(p.numel() for p in model.parameters())/1e6:.1f}")

    # Embed in batches
    embeddings = {}
    n = len(seqs)
    for i in range(0, n, args.batch):
        batch_keys = keys[i:i + args.batch]
        batch_seqs = seqs[i:i + args.batch]
        enc = tok(batch_seqs, return_tensors="pt", padding=True, truncation=True)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            out_h = model(**enc, output_hidden_states=False)
        h = out_h.last_hidden_state if hasattr(out_h, "last_hidden_state") else out_h[0]
        # Mean-pool over non-pad tokens (simple — HyenaDNA char-level doesn't need attn mask)
        emb = h.mean(dim=1).detach().cpu().to(torch.float32)
        for k, e in zip(batch_keys, emb):
            embeddings[k] = e.clone()
        if (i // args.batch) % 20 == 0:
            print(f"  {i+len(batch_keys)}/{n}")

    torch.save(embeddings, out)
    print(f"[done] wrote {out} — {len(embeddings)} embeddings, dim {next(iter(embeddings.values())).shape}")


if __name__ == "__main__":
    main()
