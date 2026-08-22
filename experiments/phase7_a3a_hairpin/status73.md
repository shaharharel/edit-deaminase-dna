
## ############################################################################
## QA TICK ~15:35 UTC — THE MIX-SHIFT CHECK. s7c prints a trinucleotide breakdown labelled
## "mix shift can imitate an editor effect" and I had never read it for this arm. It passes.
## ############################################################################

If the effect lived in ONE trinucleotide it would be a composition artefact. It does not.

    arm                ctx  stem   editor   calib   ed-cal   n_hp_ed
    Y130G-clone2       TCA     6    1.514   1.106   +0.408     1,152
                       TCA     7    2.310   1.495   +0.815       604
                       TCA     8    3.335   1.785   +1.550       309
                       TCT     6    1.209   0.996   +0.213       917
                       TCT     7    1.824   1.255   +0.569       471
                       TCT     8    2.369   1.452   +0.917       210
    Y130G-clone1       TCA     8    3.480   1.785   +1.695       292
                       TCT     8    2.775   1.452   +1.323       215
    VALIDATION (null)  TCA   6/7/8  ed-cal +0.069 / -0.075 / -0.057
                       TCT   6/7/8  ed-cal +0.031 / -0.022 / -0.048

FOUR THINGS FOLLOW
 1. THE EFFECT IS IN BOTH TRINUCLEOTIDES, both clones, every stem. A mix shift would put it
    in one context only. It is not a composition artefact.
 2. THE NULL CONTROL IS FLAT IN BOTH CONTEXTS -- ed-cal between -0.075 and +0.069 across all
    six cells. The estimator returns zero on a true-null pair in this breakdown too.
 3. TCA EXCEEDS TCT CONSISTENTLY: at stem 8, ratios 1.41 and 1.25 across the two clones. That
    ordering is what A3A biology predicts (the YTCA preference), and haA3A retains the
    APOBEC3A catalytic core.
 4. n IS AMPLE EVERYWHERE: 1,152 / 604 / 309 for TCA and 917 / 471 / 210 for TCT.

### THE CAUTION THAT MATTERS MOST HERE, AND I AM PUTTING IT IN THE SAME BREATH
Point 3 is the seductive one and I am deliberately not leaning on it. "CONSISTENT WITH A3A
BIOLOGY" IS EXACTLY THE NARRATIVE THE A3A-vs-A3B DICHOTOMY DIED FROM. That claim was
biologically plausible too, and it still went from +0.355 at n=53 to +0.079 at p=0.44 with
n=97. Structural consistency raises how INTERESTING this result is. It does not convert an
exploratory result into a finding, and it is not evidence I am entitled to spend.

### RUNNING TALLY OF CHECKS THIS POSITIVE HAS PASSED
gate A0 on both editor clones; burden matched to 1% (independently recomputed); GC-decile;
COMPLEXITY-decile; trinucleotide mix-shift; measured clone-luck floor +0.075; same-day
true-null validation flat on every endpoint and every breakdown; two-clone replication;
monotone stem scaling; adequate n at every stem, band, decile and context.
NOTHING ON THE STANDING CHECKLIST HAS MOVED IT. It is still EXPLORATORY, there is still NO
editor claim, and the inversion against A3A-Y130F is still unexplained. The result that would
change this is VA -- the OTHER haA3A variant, aligning now, ~1.2h and ~4.2h out.

### INTEGRITY -- node B clean; node A stopped deliberately with its queue complete
