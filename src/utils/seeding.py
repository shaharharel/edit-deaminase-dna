"""Deterministic seeding helper."""
from __future__ import annotations
import os, random
import numpy as np
import torch


def set_seed(seed: int = 42, deterministic: bool = True):
    """Set all seeds. If deterministic=True, also force deterministic algorithms
    (slower, especially on MPS/CUDA, but reproducible)."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        # don't enable torch.use_deterministic_algorithms by default on MPS
        # because many ops fall back to CPU at large perf cost
