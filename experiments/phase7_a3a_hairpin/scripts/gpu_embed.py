#!/usr/bin/env python
"""DeaminaFormer-DNA, stage 1: frozen embeddings from a pretrained DNA foundation model.

WHY THIS AND NOT MORE GRADIENT BOOSTING. The GB baseline sees 81 bp one-hot plus hairpin
geometry and tops out at 5.280x in the top 0.1% tail. Two facts say its ceiling is the
INPUT, not the learner: the structure branch adds +38% tail for +0.0023 AUROC, and label
purity beat label volume 4:1. A pretrained DNA LM brings 1 kb of context and representations
learned from whole genomes rather than from 84k labels.

Two encoders, chosen because they fail differently:
  NT-v2-100M     6-mer tokens, multi-species -> motif and composition grammar
  HyenaDNA-1k    single-nucleotide           -> exact spacing, palindromes, stem geometry
A 6-mer tokenizer cannot represent a 1 bp shift in a stem; for a hairpin hypothesis that is
the whole question, so both are extracted and ablated separately.

REWRITTEN after the first run died of CUDA OOM at 883,392 / 923,989 with nothing saved:

 1. AMBIGUOUS BASES BLOW UP THE TOKENIZER. NT-v2 cannot form a 6-mer across an N, so it
    falls back to single-character tokens. Measured: a clean 1,025 bp window is 176 tokens;
    N-heavy ones are 376-461; the worst window here is all N at 1,025 tokens. Padding is to
    the longest member of a batch, so ONE such window forces the whole batch to ~36x the
    attention memory -- the 6.03 GiB allocation that failed. 93 windows in 923,989 (0.0101%)
    contain any N; they are embedded one at a time, and zero-filled with an indicator if even
    that fails. Absent, never fabricated -- the same rule used for ambiguous bases elsewhere.

 2. NOTHING WAS CHECKPOINTED. Shards are written every 100k rows and merged at the end, so a
    failure costs one shard instead of the run.
"""
import os, sys, time, glob
os.environ.setdefault("HF_HOME", "/mnt/a3a/hf")
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel, AutoModelForMaskedLM

BASE = "/mnt/a3a"
WIN = f"{BASE}/a3a_trainset_v5_win1k.npz"
MODELS = {"ntv2": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
          "hyena": "LongSafari/hyenadna-tiny-1k-seqlen-hf"}
WHICH = sys.argv[1] if len(sys.argv) > 1 else "ntv2"
BATCH = int(os.environ.get("EMB_BATCH", "128"))
SHARD = int(os.environ.get("EMB_SHARD", "100000"))
LIMIT = int(os.environ.get("EMB_LIMIT", "0"))

d = np.load(WIN, allow_pickle=True)
win, chrom, pos, y = d["win1k"], d["chrom"], d["pos"], d["y"]
N = len(y) if LIMIT == 0 else min(LIMIT, len(y))
MID = win.shape[1] // 2
print(f"{WIN}: n={len(y):,}  using {N:,}  window={win.shape[1]} bp", flush=True)

LUT = np.array(list("ACGTN"))
def seqs(lo, hi):
    return ["".join(LUT[win[i]]) for i in range(lo, hi)]

bad = sum(1 for s in seqs(0, 200) if s[MID] != "C")
print(f"decode check: focal base C in {200-bad}/200 probe windows", flush=True)
if bad:
    raise SystemExit("*** decode is wrong -- focal base is not C. Nothing computed. ***")

tok = AutoTokenizer.from_pretrained(MODELS[WHICH], trust_remote_code=True)
try:
    model = AutoModel.from_pretrained(MODELS[WHICH], trust_remote_code=True); HEADLESS = True
    ENTRY = "AutoModel"
except ValueError:
    model = AutoModelForMaskedLM.from_pretrained(MODELS[WHICH], trust_remote_code=True)
    HEADLESS = False; ENTRY = "AutoModelForMaskedLM + hidden_states[-1]"
model = model.cuda().eval().half()
print(f"  entry point: {ENTRY}\n{WHICH}: {sum(p.numel() for p in model.parameters())/1e6:.0f}M params",
      flush=True)

AMBIG = (win[:N] == 4).sum(1) > 0
print(f"ambiguous windows: {int(AMBIG.sum())} / {N:,} ({100*AMBIG.mean():.4f}%) "
      f"-- run one at a time so a single bad row cannot pad a whole batch", flush=True)

SHARDDIR = f"{BASE}/shards_{WHICH}"
os.makedirs(SHARDDIR, exist_ok=True)
for f in glob.glob(f"{SHARDDIR}/*.npz"):
    os.remove(f)

buf_m, buf_c = [], []
n_zero = 0
last = 0
t0 = time.time()
with torch.no_grad():
    for lo in range(0, N, BATCH):
        hi = min(lo + BATCH, N)
        step = 1 if AMBIG[lo:hi].any() else (hi - lo)
        for a in range(lo, hi, step):
            b = min(a + step, hi)
            try:
                enc = tok(seqs(a, b), return_tensors="pt", padding=True,
                          truncation=True, max_length=1100)
                enc = {k: v.cuda() for k, v in enc.items()}
                h = (model(**enc).last_hidden_state if HEADLESS
                     else model(**enc, output_hidden_states=True).hidden_states[-1])
            except torch.OutOfMemoryError:
                torch.cuda.empty_cache()
                D = buf_m[0].shape[1] if buf_m else 512
                buf_m.append(np.zeros((b - a, D), np.float16))
                buf_c.append(np.zeros((b - a, D), np.float16))
                n_zero += b - a
                print(f"  OOM rows {a}-{b}: ZERO-FILLED, continuing", flush=True)
                continue
            am = enc.get("attention_mask")
            if am is None:
                mean = h.mean(1)
            else:
                mm = am.unsqueeze(-1).to(h.dtype)
                mean = (h * mm).sum(1) / mm.sum(1).clamp(min=1)
            buf_m.append(mean.float().cpu().numpy().astype(np.float16))
            buf_c.append(h[:, 0].float().cpu().numpy().astype(np.float16))
        if hi - last >= SHARD or hi >= N:
            np.savez(f"{SHARDDIR}/{last:09d}_{hi:09d}.npz",
                     mean=np.concatenate(buf_m), cls=np.concatenate(buf_c))
            buf_m, buf_c = [], []
            last = hi
        if lo % (BATCH * 200) == 0:
            el = time.time() - t0
            r = hi / max(el, 1e-9)
            print(f"  {hi:,}/{N:,}  {r:.0f} seq/s  ETA {(N-hi)/max(r,1e-9)/60:.0f} min", flush=True)

parts = sorted(glob.glob(f"{SHARDDIR}/*.npz"))
E_mean = np.concatenate([np.load(f)["mean"] for f in parts])
E_cls = np.concatenate([np.load(f)["cls"] for f in parts])
assert len(E_mean) == N, f"shards total {len(E_mean)} but expected {N}"
print(f"merged {len(parts)} shards, {n_zero} rows zero-filled, "
      f"{(time.time()-t0)/60:.1f} min; mean {E_mean.shape}", flush=True)
np.savez(f"{BASE}/emb_{WHICH}_v5.npz", mean=E_mean, cls=E_cls, ambiguous=AMBIG,
         chrom=chrom[:N], pos=pos[:N], y=y[:N])
print(f"wrote {BASE}/emb_{WHICH}_v5.npz", flush=True)
