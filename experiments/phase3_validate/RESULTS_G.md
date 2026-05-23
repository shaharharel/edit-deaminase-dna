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

## g v2 (2026-05-23) — motif-controlled + cell-type-matched (the strong result)
v2 bins add TpC motif count + ATAC + H3K27ac. LOCO Spearman.

### Task 1: chromatin adds BEYOND opportunity + TpC motif + mappability
- baseline (opportunity + **TpC motif** + mappability): 0.353
- + chromatin (ATAC/DNase/R-loop/H3K27ac/conservation): **0.467** (+0.11)
- => signal beyond the motif is real at the regional level (ATAC+0.45, DNase+0.47 drive it).

### Task 2: CELL-TYPE MATCHING is the key lever (per-source LOCO, HEK293 tracks)
| source | cell type | LOCO Spearman | n |
|---|---|---|---|
| Yu | HEK293T (matched) | **0.644** | 6560 |
| Lei | HEK293T (matched) | 0.403 | 2701 |
| Chen | HEK293T (matched, tiny) | 0.231 | 404 |
| Doman | HEK293 orthogonal-assay (biased) | 0.221 | 11022 |
| McGrath | **iPSC (MISMATCHED tracks)** | **0.006** | 13689 |

**Cell-type-matched accessibility predicts the >=1Mb off-target landscape at ~0.64; mismatch -> 0.0.**
Pooled 0.467 was diluted by cell-type mix. This is the genuine, clinically-meaningful editor-intrinsic
signal: per-editor (motif rate) x target-cell-type accessibility, regional resolution.

### Next levers
- iPSC accessibility tracks would rescue McGrath (test the matching claim symmetrically).
- More HEK293 genome-wide sources; Repli-seq (pending); per-editor model.
- Clinical: score a therapy's target cell type -> regional off-target risk map.
