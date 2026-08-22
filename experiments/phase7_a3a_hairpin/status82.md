
## ############################################################################
## QA TICK ~19:55 UTC — CHECK 1 ON THE LEARNABILITY TEST. I had ASSERTED the matching
## "by construction" instead of measuring it. Measured now. It holds exactly.
## ############################################################################

### WHY I DID NOT TRUST MY OWN ASSERTION
s13_learnable.py printed "negatives matched on trinuc AND strand BY CONSTRUCTION". That phrase
is what I also said about the atomic write (which silently mis-named 23 files) and about the
gate's peer list (which silently fell back to 2 peers on node B). It is not evidence.
AND THERE WAS A SPECIFIC WAY IT COULD BE FALSE. The matcher does
    need = (#positives in this tri x strand stratum) * 10
    take = min(need, len(candidates))
so if ANY stratum ran short of candidates, the ratio silently drops for that stratum alone,
the class balance breaks, and NOTHING WARNS.

### MEASURED, BOTH ARMS, ALL FOUR STRATA
    editor        tri  strand  positives  negatives  ratio   pos share  neg share
                  TCA       0     21,390    213,900  10.000     0.2561     0.2561
                  TCA       1     21,553    215,530  10.000     0.2580     0.2580
                  TCT       0     20,275    202,750  10.000     0.2427     0.2427
                  TCT       1     20,312    203,120  10.000     0.2432     0.2432
    calibrator    TCA       0     21,575    215,750  10.000     0.2560     0.2560
                  TCA       1     21,817    218,170  10.000     0.2589     0.2589
                  TCT       0     20,487    204,870  10.000     0.2431     0.2431
                  TCT       1     20,401    204,010  10.000     0.2421     0.2421
    MAX |pos share - neg share| = 0.00000 in BOTH arms
    strata that ran short of candidates: 0 in BOTH arms -- the 10:1 ratio held everywhere

### SO CHECK 1 PASSES ON THE LEARNABILITY TEST
Focal base is C and offset -1 is T for 100% of both classes by the universe definition (it is
a TCW universe). Offset +1 is the W, which is exactly the `tri` field, and its distribution is
IDENTICAL between classes to five decimal places in all four strata. Strand likewise. The
4.422x / 2.350x numbers rest on a genuinely matched design, not an asserted one.

### THE POINT OF DOING THIS
The learnability result is the one that answered the user's question, and the one on which I
based the statement that the obstacle is the DATA rather than the model. That statement is only
as good as the matching underneath it. It would have been cheap to leave "by construction" in
place; it cost one script to know instead.

### INTEGRITY -- node B clean, node A stopped with its queue complete
