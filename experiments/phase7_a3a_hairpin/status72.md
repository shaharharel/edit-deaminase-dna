
## ############################################################################
## QA TICK ~14:55 UTC — CHECK 4's UNPERFORMED HALF. Complexity stratification had NEVER
## been run on anything. It is the check most likely to bite this signal. It passes.
## ############################################################################

The standing checklist requires signals to survive GC-decile AND COMPLEXITY stratification.
GC has been done exhaustively all session; COMPLEXITY had never been done at all -- and it is
the check most likely to kill THIS signal, because hairpins ARE inverted repeats and inverted
repeats live in low-complexity sequence where mappability and alignment are worst.

### THE MEASURE, AND WHY IT IS NOT CIRCULAR
No raw sequence is stored in the universe files, so k-mer entropy is unavailable. I used LOCAL
TCW SITE DENSITY -- universe sites within +/-500bp -- which packs tightly in repeats and
spreads in ordinary sequence. It is computed from `pos` alone and is INDEPENDENT of every
hairpin feature, so unlike stratifying on stem length it is not circular.

### RESULT: THE SIGNAL IS UNIFORM ACROSS COMPLEXITY, NOT CONCENTRATED IN REPEATS
    dec  editor   calib   ed-cal   n_hp
      0   1.189   1.088   +0.101    176   <- least repetitive
      1   1.245   1.039   +0.206    185
      2   1.263   0.979   +0.284    155
      3   1.529   1.001   +0.528    216
      4   1.333   1.074   +0.259    201
      5   1.316   1.013   +0.303    194
      6   1.397   0.976   +0.421    188
      7   1.473   1.283   +0.190    166
      8   1.304   0.846   +0.458    172
      9   1.468   0.953   +0.515    230   <- most repetitive
    MANTEL-HAENSZEL:  editor 1.358   calibrator 1.022   ed-cal +0.336
    GC-stratified:    editor 1.357   calibrator 1.039   ed-cal +0.318   ESSENTIALLY IDENTICAL
The editor is above the calibrator in 10/10 deciles and clears the +0.11 bar in 9/10, with
n_hp 155-230 in every stratum -- no thin cells. The calibrator sits at ~1.0 throughout. If
this were an alignment artefact of low-complexity inverted repeats it would CONCENTRATE in the
high-density deciles. It does not; it is present everywhere.

### AND THE CAVEAT I AM NOT POOLING AWAY
There IS a complexity trend: corr(decile, ed-cal) = +0.594, and DECILE 0 -- the least
repetitive, most complex sequence -- gives ed-cal +0.101, INSIDE the +0.11 bar and near the
+0.075 clone-luck floor. The effect is weakest exactly where complexity is highest. It does
not vanish there, and it does not concentrate in the repeat tail either, but the gradient is
real and a reader should see it rather than only the pooled +0.336.

### STATUS OF THE POSITIVE AFTER THIS
Checks now passed: gate A0 on both editor clones, burden matched to 1%, GC-decile, COMPLEXITY-
decile, measured clone-luck floor, same-day true-null validation, two-clone replication,
monotone stem scaling, and adequate n at every stem and every stratum.
It is still EXPLORATORY under its pre-registration and I am still making NO editor claim. The
inversion against A3A-Y130F is still unexplained. What I can say is that I have now run every
check on the standing list against it and none of them has moved it.

### HEALTH
ai-chem  TERMINATED by my own stop, queue complete. ai-chem2 3 bwa, 0 failed, 196G,
VA-clone1 69.2% (1.7h), VA-clone2 25.2% (4.7h), YE1 x2 queued.
