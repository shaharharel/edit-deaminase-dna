"""Splitters for DeaminaFormer training.

QA I16 fix: random per-site splits leak batch effects across PRJNAs.
ONLY leave-one-PRJNA-out is safe for Levanon experiments.
"""
from __future__ import annotations
from typing import Iterator
import numpy as np
import pandas as pd


def leave_one_prjna_out(df: pd.DataFrame, prjna_col: str = "dataset"
                          ) -> Iterator[tuple[np.ndarray, np.ndarray, str]]:
    """Yield (train_idx, test_idx, held_out_prjna) tuples for each PRJNA in df.

    Use this — NOT random splits — for all Levanon training/eval.
    """
    prjnas = sorted(df[prjna_col].unique())
    for prjna in prjnas:
        test_mask = (df[prjna_col] == prjna).values
        train_idx = np.where(~test_mask)[0]
        test_idx = np.where(test_mask)[0]
        yield train_idx, test_idx, prjna


def leave_one_editor_out(df: pd.DataFrame, editor_col: str = "editor"
                          ) -> Iterator[tuple[np.ndarray, np.ndarray, str]]:
    """Yield (train_idx, test_idx, held_out_editor) — tests enzyme-token generalization.
    Use SPARINGLY (we have few editors; QA I12 says this is mostly diagnostic, not validation).
    """
    editors = sorted(df[editor_col].unique())
    for editor in editors:
        test_mask = (df[editor_col] == editor).values
        train_idx = np.where(~test_mask)[0]
        test_idx = np.where(test_mask)[0]
        yield train_idx, test_idx, editor


def stratified_negative_sample(df: pd.DataFrame, n_negatives_per_positive: int = 1,
                                  match_columns: list[str] | None = None,
                                  hotspot_col: str = "hotspot",
                                  seed: int = 42) -> pd.DataFrame:
    """Sample matched negatives for each positive (hotspot=1) row.

    QA I17 fix: ensure BE samples have BOTH positive cytidines AND negative cytidines
    so the model can't shortcut via substrate token.

    match_columns: list of columns to match on (e.g., ["editor", "coverage_decile"]).
    """
    rng = np.random.default_rng(seed)
    if match_columns is None:
        match_columns = ["editor"]
    pos = df[df[hotspot_col] == 1]
    neg_pool = df[df[hotspot_col] == 0]
    out_negs = []
    for _, p_row in pos.iterrows():
        cands = neg_pool
        for col in match_columns:
            cands = cands[cands[col] == p_row[col]]
        if len(cands) == 0:
            continue
        n = min(n_negatives_per_positive, len(cands))
        out_negs.append(cands.sample(n=n, random_state=rng.integers(0, 1_000_000)))
    if not out_negs:
        return pd.concat([pos, pd.DataFrame(columns=df.columns)], ignore_index=True)
    sampled_negs = pd.concat(out_negs, ignore_index=True)
    return pd.concat([pos, sampled_negs], ignore_index=True)
