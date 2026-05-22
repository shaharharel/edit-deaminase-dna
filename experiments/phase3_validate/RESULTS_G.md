# g-landscape — does chromatin predict the >=1Mb off-target landscape? (2026-05-23)

Per-1Mb bins (2,838), guide-indep CBE edit density vs chromatin tracks. Leave-one-chromosome-out.

## Raw Spearman with edit_count
DNase **+0.47**, R-loop **+0.39**, eligible-C +0.38, mappability +0.13, phyloP +0.11, phastCons +0.12.

## Label reproducibility at 1Mb (target ceiling)
Yu↔Lei +0.26, Lei↔Chen +0.28, Yu↔Chen +0.15 (genome-wide assays agree); Doman ~0 (orthogonal R-loop
assay = guide-placed), McGrath ~0 (clonal). => achievable correlation is bounded by ~0.26 label noise.

## Leave-one-chromosome-out Spearman(pred, observed)
| model | held-out Spearman |
|---|---|
| baseline (eligible-C + mappability) | 0.351 |
| + R-loop | 0.375 |
| + chromatin (R-loop, DNase, conservation) | **0.419** |

## Conclusion
**Chromatin predicts the megabase guide-independent off-target landscape beyond opportunity+ascertainment**
(+0.07 LOCO, generalizes across chromosomes), driven by DNase (open chromatin) + R-loop (ssDNA exposure).
This is the genuine ML signal (non-motif, learnable) — unlike the trivial f preference.
MODEST effect, bounded by low label reproducibility. NEXT: add replication timing (Repli-seq, the key
ssDNA driver, not yet sourced) + fix histone tracks (.bigWig); expect further lift.
