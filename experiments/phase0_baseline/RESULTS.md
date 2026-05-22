# P0 — Deaminase-preference transferability (2026-05-22)

Sequence-only (±10 bp local context, center dropped), numpy logistic regression.
Positives = oriented edited-C contexts; negatives = random genomic C (1:1).
Question: does guide-independent deaminase preference learned on one source/cell-type
transfer to another?

## Transfer AUROC (rows=train, cols=test)

| train\test | Doman | McGrath | Yu | Lei | Chen |
|---|---|---|---|---|---|
| Doman (rAPOBEC1, HEK293) | 0.824 | **0.820** | 0.521 | 0.482 | 0.162 |
| McGrath (rAPOBEC1, iPSC) | **0.798** | 0.839 | 0.531 | 0.539 | 0.381 |
| Yu (eng. YE1) | 0.537 | 0.514 | 0.855 | 0.115 | 0.209 |
| Lei (BE4max?) | 0.518 | 0.525 | 0.352 | 0.989 | 0.807 |
| Chen (dual tdCBE) | 0.322 | 0.348 | 0.425 | 0.862 | 0.985 |

Within-source 5-fold ceiling: Doman 0.819, McGrath 0.832, Yu 0.848, Lei 0.987(!), Chen 0.968(!).

## Conclusions
1. **rAPOBEC1 preference transfers across cell type/lab**: Doman↔McGrath 0.82/0.80 ≈ ceilings.
   The transferable f_deaminase factor is real and site-learnable. HYPOTHESIS CONFIRMED.
2. **Preference is editor-specific**: rAPOBEC1 model fails on Yu (eng.), Chen (dual) → the
   decomposed deaminase embedding is required.
3. **Lei is a confirmed data artifact** (3 independent signals): CCC motif (not TC),
   no transfer from genuine rAPOBEC1 (0.48), inflated within-CV (0.987). FIX or DROP before training.
4. Sequence-only ceiling ~0.83 for clean rAPOBEC1 = the f factor; accessibility (g)
   features add the cell-type-specific component for full burden.
