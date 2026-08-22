
## ############################################################################
## 2026-08-22 19:27 UTC — VA-clone2 LANDED POSITIVE AND DISAGREES WITH VA-clone1 IN SIGN.
## A FOURTH OUTCOME I DID NOT PRE-REGISTER, AND IT DAMAGES THE Y130G RESULT.
## ############################################################################

    arm                     GC-adj MH ed      cal    ed-cal   complexity ed-cal   category
    VA-clone1 vs nCas9-c1          0.975    1.041    -0.066              -0.066   NULL
    VA-clone2 vs nCas9-c1          1.215    1.042    +0.173              +0.197   POSITIVE
    Y130G-clone2 vs nCas9          1.357    1.039    +0.318              +0.336   POSITIVE
    Y130G-clone1 vs nCas9          1.437    1.041    +0.396                  --   POSITIVE
VA-clone2 banded MH: stem6 1.215/1.064, stem7 1.682/1.365, stem8 2.311/1.629.

### THE PRE-REGISTRATION DID NOT ANTICIPATE THIS
At 15:35 I wrote three outcomes for VA -- positive means class effect, null means
Y130G-specific, negative means the endpoint is broken -- and committed to reporting whichever
landed without reinterpreting. VA-clone1 is the NULL category. VA-clone2 is the POSITIVE
category. A FOURTH OUTCOME LANDED: THE CLONES DISAGREE.
Neither category cleanly applies. I am not picking the one I prefer, and I am not retrofitting
a fourth category that happens to be convenient. What I can say is what was measured.

### AND THIS IS THE PART THAT HURTS
    same-editor clone spread, Y130G:  0.078   (+0.318 vs +0.396)
    same-editor clone spread, VA:     0.239   (-0.066 vs +0.173)
    the clone-luck floor I have quoted ALL SESSION: +0.075
THAT FLOOR WAS MEASURED ON Y130G-clone1 vs Y130G-clone2 -- THE TIGHTER OF THE TWO PAIRS.
VA's same-editor spread is 3.1x Y130G's and 3.2x the floor I have been comparing everything
against. If 0.239 is the true same-editor variability, Y130G's mean of +0.357 sits 1.5 SIGMA
above it, not the 4.8 sigma it looked like against +0.075.

### I CANNOT DISMISS THIS AS "VA IS JUST NOISY"
Both VA clones PASSED gate A0 on the corrected depth-adjusted criterion. Both ran against the
same calibrator, through the same audit chain, on the same day. There is no quality basis for
discarding either one, and discarding the inconvenient one is exactly the move this project's
QA discipline exists to prevent.

### WHAT THIS DOES TO THE RECORD
Sections 25 / 25.1-25.5 reported the Y130G positive as having survived every check, with the
clone-luck floor of +0.075 as the yardstick. THAT YARDSTICK IS NOW IN DOUBT, because a second
variant's clones -- equally qualified -- differ by three times as much. The Y130G effect is
still the largest measured and still replicated across its own two clones, but the margin over
same-editor variability is much smaller than I have been reporting.
THE HONEST POSITION IS NOW WEAKER THAN IT WAS AN HOUR AGO, and it moved because more data
arrived, which is the only reason a position should move.

### WHAT WOULD SETTLE IT
More clones per variant. Two is not enough to estimate same-editor variability -- that is what
this result demonstrates, and it is the same bound the clonal-power check gave from the other
direction (25-50 clones per arm). YE1 x2 is aligning and will add a third variant, but it will
add another PAIR, not the spread estimate that is actually needed.
