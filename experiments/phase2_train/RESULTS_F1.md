# F1 — Full model vs simple probe (2026-05-23)

Task: held-out transfer (train Doman → test McGrath, rAPOBEC1), RANDOM-negative regime
(0.82-comparable). Feature ablation, MLP, MPS. HyenaDNA-small embedded on V100 (idle), 26,071 seqs.

| Feature set | dim | Doman→McGr | McGr→Doman |
|---|---|---|---|
| seq ±10bp one-hot | 80 | **0.819** | **0.792** |
| handcrafted-seq (trinuc+GC+kmer) | 83 | 0.817 | 0.769 |
| HyenaDNA mean-pool (513bp) | 256 | 0.560 | 0.555 |
| HyenaDNA center-token (±16) | 256 | 0.638 | 0.629 |
| seq10 + HyenaDNA-center | 336 | 0.802 | 0.773 |
| all | 419 | 0.805 | 0.773 |

## Conclusion
**The simple local-sequence model wins decisively.** Frozen HyenaDNA embeddings (mean OR
center-pooled) underperform ±10bp one-hot by ~0.18 AUROC and do not help when concatenated.
Mechanistic reason: guide-independent deaminase preference is a sharp ±2bp local motif (TpC);
a long-range genome-LM embedding is a worse representation of local chemistry.

## Implication for Evo-2 / foundation models
For the **f (deaminase chemistry)** factor, foundation models add no value (signal is local).
Evo-2 (larger, longer-context) would likely also not beat one-hot for f. Foundation models are
only plausibly useful for the **g (long-range chromatin/accessibility)** factor — a different,
harder problem with its own validity issues (location reproducible only ≥1 Mb). Recommendation:
keep the simple sequence model for f; do NOT invest in FM embeddings for the preference model.
