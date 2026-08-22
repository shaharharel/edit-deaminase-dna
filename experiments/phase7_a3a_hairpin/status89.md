
## ############################################################################
## 2026-08-22 ~21:40 UTC — THREE-PANEL SYNTHESIS. The decisive test they ranked #1 CLEARS
## Y130G; three other findings damage it more; and the literature reframes the whole null.
## ############################################################################

### 1. THE SHARED-ANCESTOR TEST -- RANKED #1 BY THE PANEL -- REFUTES ITS OWN HYPOTHESIS
The cross-clone recurrence filter had NEVER been applied to Y130G. It killed the last two
positives in this project (A3A-Y130F stem-8 went 1.326 shared -> 1.026 private). Applied now,
removing 80,366 sites from ELIGIBILITY (both numerator and background, as the patch intends):
                    stem6         stem7         stem8      ed-cal
    WITHOUT     1.411/1.086   2.101/1.396   2.969/1.660   +0.325 +0.705 +1.309
    WITH        1.407/1.076   2.078/1.347   2.945/1.587   +0.331 +0.731 +1.358
THE EFFECT IS UNCHANGED, MARGINALLY LARGER. Y130G's signal is NOT a shared pre-isolation
subclonal ancestor. The panel's strongest alternative explanation is refuted by its own test.

### 2. BUT THE FEASIBILITY CONTRADICTION STANDS (section 38, verified)
230 excess stem>=8 calls require 12-243x hairpin selectivity at every M the burden endpoint
permits. Endogenous A3A tops out at 2-3x. Empty intersection.

### 3. MY STATISTICS WERE WRONG IN A WAY THAT CHANGES THE ANSWER
Label permutation over the 6 within-study clones: P(both Y130G clones rank top-2) = 1/15 =
0.067; x3 for choosing Y130G post hoc = ~0.20. p ~ 0.07 UNCORRECTED, ~0.2 CORRECTED.
SUGGESTIVE, NOT SIGNIFICANT. That replaces the entire "1.5 to 4.8 sigma" discussion, which was
never licensed -- both ends are |difference| of two clones, df=1, t_0.975 = 12.7.

### 4. I HAD FIVE WITHIN-GENOTYPE SPREADS IN MY OWN RECORD, NOT TWO
    A3A-Y130F c2 vs c5   0.021   editor
    nCas9 c1 vs c2       0.036   DEAMINASE-FREE
    Y130G c1 vs c2       0.075   editor
    D10A c1 vs c6        0.117   DEAMINASE-FREE   <-- EXCEEDS MY +0.11 PRE-REGISTERED BAR
    VA c1 vs c2          0.237   editor
An 11x span, and A NO-DEAMINASE CLONE PAIR CLEARS THE BAR. A bar a placebo pair exceeds is not
a bar. Quoting a single "clone-luck floor" was a modelling choice, not a measurement.

### 5. THE LITERATURE REFRAMES THE NULL ENTIRELY
DOMAN 2020 RAN THIS EXACT EXPERIMENT AT n=8/8/7 (BE4/YE1/nCas9) AND YE1 WAS A WGS NULL.
YE1's 10-100x superiority was established by the ORTHOGONAL R-LOOP ASSAY, not by WGS. Our null
is the literature-expected result: we reproduced Doman's YE1 arm across four engineered
variants. Doman himself states stochastic Cas9-independent deamination is "well below the
~0.1% detection limit of practical high-throughput DNA sequencing".

AND HEK293T IS CLOSE TO THE WORST POSSIBLE HOST:
  - it expresses endogenous APOBEC3B (A3A not detectable by immunoblot)
  - PLASMID TRANSFECTION ITSELF INDUCES A3B, published, rising 24-72h -- the background
    deaminase is switched on by the act of delivering the editor
  - APOBEC mutagenesis in cell lines is EPISODIC BETWEEN SUBCLONES (Petljak, Cell 2019 /
    Nature 2022) -- a mechanism for VA's sign flip that needs no editor at all
AND: NO PUBLISHED METHOD SEPARATES EDITOR-DEPOSITED FROM ENDOGENOUS APOBEC BY SIGNATURE ALONE.
Our editors are A3A-DERIVED. We have been trying to distinguish A3A-family editors from
endogenous A3A/A3B using the A3A signature. The axis is a priori non-discriminative.

### 6. THE NEXT ACTION, AND IT IS CHEAP
RUN DOMAN'S BE4 (PRJNA553240, n=8, KNOWN-POSITIVE ARM) THROUGH OUR EXACT PIPELINE. Same cell
line, same assay, same clone lottery. IT IS THE POSITIVE CONTROL THIS ENTIRE SESSION HAS
LACKED. If the pipeline recovers BE4 > nCas9 and still fails on YE1/Y130G, the null is real and
bounded. If it cannot recover BE4, we have a pipeline problem, not a biology result.
Everything else the panel ranked -- more clones (565-9,000/arm), regional aggregation, bigger
models, DNA LMs -- they rated not worth doing, and the arithmetic behind each is in sections
33-38.
