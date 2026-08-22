
## ############################################################################
## 2026-08-22 ~17:05 UTC — THE VA PANEL TEST CLOSES THE MECHANISM QUESTION.
## The transfer is CONDITIONAL, and the condition was predicted before the test ran.
## ############################################################################

### VA-clone1 panel, PCAWG-trained model, 200-permutation null
    panel   n_ed    ed x  n_cal   cal x      gap   null sd      z
    0.01%      9   0.954     15   1.774   -0.820     0.470   -1.8
    0.10%    113   1.198    114   1.348   -0.151     0.155   -1.0
    1.00%   1018   1.079    944   1.117   -0.038     0.047   -0.8
    5.00%   4618   0.979   4077   0.964   +0.014     0.020   +0.7
NULL AT EVERY DEPTH, |z| < 2 throughout. ed-specific 94,361 vs cal-specific 84,542, so 1.12x
matched -- not a burden artefact.

### THE THREE-ARM PICTURE IS NOW COMPLETE AND INTERNALLY CONSISTENT
    editor       hairpin excess         panel transfer
    Y130G        +0.318 / +0.396        +0.466 / +0.459   z=2.7  (z=5.0 at 0.01%)
    VA           -0.066                 -0.151            z=-1.0   NULL
    A3A-Y130F    -0.18  / -0.21         -0.483            anti-transfer
THE MODEL TRANSFERS IF AND ONLY IF THE TARGET ARM CARRIES THE PROPERTY THE MODEL LEARNED.
The transfer is not sporadic; it is CONDITIONAL on hairpin excess being present, and the
prediction was recorded BEFORE the confirming test was run.

### THIS ALSO DISSOLVES THE "BIOLOGICAL INVERSION" THAT HAS BEEN OPEN SINCE SECTION 25
Sections 25 / 25.1 / 26 / 27 / 28 all circled the same puzzle: the engineered-safe variant
showed signal and the active one did not, which looked backwards. IT WAS NEVER AN INVERSION IN
THE ENDPOINT. Y130G genuinely carries hairpin-structured off-targets; VA and A3A-Y130F do not;
and the PCAWG model tracks that faithfully in all three directions INCLUDING THE NEGATIVE ONE.
An endpoint that produced a positive, a null and an anti-transfer in exactly the three arms
where hairpin excess is present, absent and absent is not a broken endpoint.

### WHAT IT STILL DOES NOT BECOME
Burden prediction. The z=5.0 panel captures 35 of 83,530 editor mutations -- it misses 99.96%
-- and the skill is gone by the panel size where capture becomes meaningful. Enrichment and
capture remain in different places, and that gap is untouched by anything above.

### THE DEFENSIBLE CLAIM, IN ONE SENTENCE
A cancer-trained hairpin model ranks base-editor off-target sites in the one editor variant
whose off-targets are hairpin-structured, is null in a variant that is not, and anti-transfers
in a third -- a conditional, mechanistic transfer, demonstrated in three arms with the
condition stated in advance.
