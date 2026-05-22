"""DeaminaFormer-DNA: BE DNA off-target safety gate model.

Architecture B (per expert advisor report):
  - Sequence encoder: Evo-2 1B embeddings (1024-d, pre-computed on A100)
  - Structural features: ENCODE per-site (DNase, R-loop, ATAC, mappability, region class)
  - Shared enzyme/BE-config embedding (same CATALOG as DeaminaFormer-RNA)
  - Substrate token: DNA_R_loop or DNA_replication_fork
  - Binary classification head (off-target probability) + optional editor-aware head

Compared to DeaminaFormer-RNA:
  - swap RNA-FM 640d -> Evo-2 1024d
  - swap loop 9d -> ENCODE features (~6d to start)
  - same hybrid editor embedding (locked EDITOR_ID_TO_IDX from src.data.editor_features)
"""
from __future__ import annotations
import sys
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

# Reuse the shared editor catalog from edit-deaminase repo
sys.path.insert(0, "/Users/shaharharel/Documents/github/edit-deaminase")
from src.data.editor_features import CATALOG, EditorEmbedding  # type: ignore

# DNA substrate tokens (subset of edit-deaminase's full SUBSTRATES list)
DNA_SUBSTRATES = [
    "DNA_R_loop",            # active R-loop window (DRIP-seq +)
    "DNA_replication_fork",  # near replication origin / TADs
    "DNA_open_chromatin",    # ATAC-positive
    "DNA_quiescent",         # ATAC-negative, low expression
]
N_DNA_SUBSTRATES = len(DNA_SUBSTRATES)

D_EVO2 = 1024              # Evo-2 1B embedding dim (TBD — confirm at first compute)
D_ENCODE = 6               # dnase, rloop, atac_overlap, mappability, gc, region_class_id
D_EDITOR = 128
D_HIDDEN = 256


class DeaminaFormerDNA(nn.Module):
    """Multi-modal DNA off-target predictor.

    Forward inputs:
      seq_emb: [B, D_EVO2] — Evo-2 mean-pooled embedding around site
      encode_feats: [B, D_ENCODE] — DNase / R-loop / ATAC / mappability / GC / region
      editor_cfgs: List[EditorConfig] — length B, looked up in shared CATALOG
      substrate_ids: [B] long tensor — index into DNA_SUBSTRATES

    Output:
      logit: [B] — binary off-target probability logit
      z: [B, D_HIDDEN] — final hidden state (for downstream calibration / ablation)
    """

    def __init__(self, d_evo2: int = D_EVO2, d_encode: int = D_ENCODE,
                  d_editor: int = D_EDITOR, d_hidden: int = D_HIDDEN, dropout: float = 0.3):
        super().__init__()
        self.seq_proj = nn.Sequential(
            nn.Linear(d_evo2, d_hidden), nn.GELU(), nn.Dropout(dropout))
        self.enc_proj = nn.Sequential(
            nn.Linear(d_encode, d_hidden), nn.GELU(), nn.Dropout(dropout))
        self.substrate_emb = nn.Embedding(N_DNA_SUBSTRATES, d_hidden)
        self.editor_emb = EditorEmbedding(d_out=d_editor)
        # Combine 4 streams via small cross-attention with a learned CLS query
        self.cls = nn.Parameter(torch.randn(1, 1, d_hidden))
        self.attn = nn.MultiheadAttention(d_hidden, num_heads=4, batch_first=True,
                                          dropout=dropout)
        self.editor_lift = nn.Linear(d_editor, d_hidden)
        self.head = nn.Sequential(
            nn.LayerNorm(d_hidden), nn.Linear(d_hidden, d_hidden), nn.GELU(),
            nn.Dropout(dropout), nn.Linear(d_hidden, 1))

    def forward(self, seq_emb, encode_feats, editor_cfgs, substrate_ids):
        # Project each modality to d_hidden
        s = self.seq_proj(seq_emb)                  # [B, H]
        e = self.enc_proj(encode_feats)             # [B, H]
        sub = self.substrate_emb(substrate_ids)      # [B, H]
        ed = self.editor_lift(self.editor_emb(editor_cfgs))  # [B, H]

        # Stack as a [B, 4, H] sequence and attend with CLS
        tokens = torch.stack([s, e, sub, ed], dim=1)  # [B, 4, H]
        B = tokens.shape[0]
        cls = self.cls.expand(B, -1, -1)              # [B, 1, H]
        out, _ = self.attn(cls, tokens, tokens)
        z = out.squeeze(1)                            # [B, H]
        logit = self.head(z).squeeze(-1)              # [B]
        return {"logit": logit, "z": z}
