# Showcase plan — confound-proof demonstration of guide-independent off-target prediction (cohort design)

Gene density (r=0.61 with landscape) is the adversary; it caused 3 of 4 over-claims. The showcase must
hold gene density constant BY DESIGN, not by regression.

## Task hierarchy
1. PRIMARY (control-by-design, no new data): within-gene-density-stratum prediction. DONE: pooled 0.351,
   DNase-alone ~0.25 within strata -> accessibility predicts at constant gene density.
2. HEADLINE (needs paired WGS): predict the TREATED-MINUS-CONTROL excess landscape. Subtraction cancels
   gene-density/mappability/germline ascertainment (shared by treated+control) -> remaining = editor-attributable.
3. SPECIFICITY (needs >=4-5 distinct-motif editors): per-editor predicted x observed MATRIX; diagonal must
   beat off-diagonal (retires the circular cross-class claim).
4. TRUMP CARD (wet lab): prospective targeted deep-seq of gene-density-MATCHED high- vs low-rank panels.

## New data to acquire (priority)
- Paired treated+CONTROL bulk WGS (non-negotiable for excess landscape).
- >=4-5 editors spanning distinct motifs (BE4max, YE1, A3A-BE, eA3A, ABE8e, hi-fi ABE).
- >=2 divergent cell types (HEK293 anchor + HSPC/hepatocyte) for the falsification test.
- >=30-40x depth both arms; replicates for per-editor noise ceiling.
- AVOID: clonal/de-novo (dispersed, unpredictable) and orthogonal/guide-placed (artifact) for the landscape.

## Honest framing
Opportunity+gene-density dominate raw 0.63 (NOT the contribution). The contribution = isolating & validating
the accessibility increment (within-stratum 0.35 / DNase 0.25; +0.04 partial). Report within-stratum + excess
numbers, not raw 0.63. Negatives are findings (cell-type-invariance => no cell-type specificity at 1Mb;
clonal assays unpredictable). Defensible claim: a >=1Mb relative-rank screening prior to triage regions.
