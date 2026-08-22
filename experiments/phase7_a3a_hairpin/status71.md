
## ############################################################################
## QA TICK ~14:15 UTC — I SUSPECTED A DOUBLE STANDARD IN MY OWN REPORTING AND CHECKED IT.
## It was not one. But the check surfaced band-level texture I should have shown.
## ############################################################################

I refused to quote the stem-8 number for the A3A-Y130F arm on the grounds that its banded
cells were too thin, and then quoted stem 8 for Y130G as the LARGEST effect (+1.309/+1.570)
WITHOUT ever checking its n. That looked like applying a rule only when it suited the
conclusion. Checking it:
    arm                  stem6 n_hp   stem7 n_hp   stem8 n_hp
    A3A-Y130F-clone2            302          119           43   <- stem 8 NOT quoted
    Y130G-clone2              2,069        1,075          519   <- stem 8 quoted
Y130G's stem 8 carries 12x the hairpin hits of the arm whose stem 8 I declined to quote. THE
STANDARD WAS CONSISTENT. 43 is thin; 519 is not.
Supporting, from the same output with p_ed = 0.0000 at every stem:
    stem 6  observed 1.362  null95 0.95-1.05   outside by +0.312
    stem 7  observed 2.069  null95 0.91-1.09   outside by +0.979
    stem 8  observed 2.866  null95 0.86-1.14   outside by +1.726

### BUT THE CHECK FOUND SOMETHING I SHOULD HAVE REPORTED AND DID NOT
The pooled MH hides the band structure. At stem 8:
    band            editor   calib   n_hp ed   n_hp cal
    (8,15)           3.269   1.955         8          5   thin
    (15,25)          3.681   1.571       110         50
    (25,35)          2.693   1.831        70         51
    (35,60)          2.503   1.600        31         21
    (60,100)         0.000   1.686         0          2   ZERO CELL
The effect lives in (15,25) and (25,35), which is where MH puts its weight. The deep bands are
thin, one is a ZERO CELL, and AT STEM 7 THE DEEPEST BAND INVERTS (editor 0.907 vs calibrator
1.339). None of that changes the pooled number -- MH weights by cell size -- but quoting only
the pooled value concealed it, and a reader is entitled to the texture. Recording it now.

### WHAT I AM STILL NOT SAYING
The positive remains EXPLORATORY under its pre-registration. Nothing here converts it. The
band structure is consistent with a real effect concentrated where the data are (mid-depth
bands) and noisy where they are thin -- but "consistent with" is not "demonstrates", and the
inversion against A3A-Y130F is still unexplained.

### INTEGRITY -- node B clean; node A stopped (deliberately, queue complete)
