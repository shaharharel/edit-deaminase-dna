
## ############################################################################
## QA TICK ~18:45 UTC — ATTACKED SECTION 34 AND FOUND A REAL LEAK. Fixed the generator,
## reran, AND THE FIX CORRECTS MY OWN OVER-STATEMENT: adjusted, not invalidated.
## ############################################################################

### THE ATTACK, AND IT WAS WARRANTED
s14 fitted ONE model on all of v5 and scored the genome, then joined mutations from THE SAME
21 DONORS the model trained on. Measured the overlap exactly, by (chrom, pos):
    mutations landing on the TCW universe (s14's denominator)   92,823
    of those, ALSO v5 training positives                        83,977
    IN-SAMPLE FRACTION                                          90.5%
And the comparison to the editor was structurally unfair: the editor arm was scored by a model
trained on PCAWG that had never seen a HEK293T site -- genuinely out-of-sample -- while the
cancer arm was scoring its own training data.

### FIXED AT THE GENERATOR AND RERUN
s14b: held-out CHROMOSOME folds, the rule every other model here obeys. Train on 4/5 of
chromosomes, score the universe sites of the held-out fifth, assemble genome-wide OUT-OF-FOLD
scores. The model that scores a site never trained on its chromosome.
                     top0.01%  top0.1%   top1%   top5%  top10%
    in-sample  (v1)    20.146    7.337   3.169   1.796   1.477
    OUT-OF-FOLD (v2)   16.698    6.819   3.007   1.665   1.372
    RANDOM     (v2)     1.077    1.056   1.044   1.001   1.003

### I SAID THE LEAK "INVALIDATES SECTION 34". THAT WAS TOO STRONG AND I AM CORRECTING IT.
The leak inflated the numbers by 7-17%, not fundamentally. Raising the concern was right and
fixing the generator was right; the measured impact is smaller than I asserted while the fix
was still running. I should not have used "invalidates" before I had the corrected number.
WHY SO LITTLE INFLATION AT 90.5% OVERLAP: the model has SIX STRUCTURAL FEATURES and no site
identity. There is nothing to memorise. It learns a hairpin-geometry function and that
function generalises across chromosomes -- which is itself evidence the signal is geometry
rather than remembered positions.

### THE CORRECTED LIKE-FOR-LIKE COMPARISON, BOTH ARMS NOW OUT-OF-SAMPLE
    panel   CANCER-OOF cap  CANCER-OOF enr   EDITOR cap  EDITOR enr   ratio
    0.10%        0.68%          6.819           0.17%      1.748      3.90x
    1.00%        3.01%          3.007           1.15%      1.146      2.62x
    5.00%        8.33%          1.665           4.93%      0.985      1.69x
Base rates within 2% (1 in 2,548 vs 1 in 2,574). Cancer concentrates 3.90x better than the
editor at top 0.1% and captures 13.72% of APOBEC-context mutations in the top 10% of sites,
while the editor's signal decays to 0.985x by top 5%.

### WHAT STANDS
The conclusion is unchanged and now rests on a fair test: THE BIOLOGY IS FINDABLE IN A WHOLE
GENOME, THE EDITOR'S 0.5-2% SLICE OF IT IS NOT. The numbers behind it are 6.819x / 3.007x /
1.665x, not the 7.337x / 3.169x / 1.796x I reported an hour ago.
