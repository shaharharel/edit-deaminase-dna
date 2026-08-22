
## ############################################################################
## 2026-08-22 ~21:00 UTC — EXPERT PANEL FOUND A DECISIVE ARGUMENT I MISSED.
## THE Y130G POSITIVE IS ARITHMETICALLY INCOMPATIBLE WITH ITS OWN BURDEN ENDPOINT.
## Sections 25 / 25.1-25.5 must be read with this attached.
## ############################################################################

### THE FEASIBILITY CONTRADICTION, VERIFIED INDEPENDENTLY
    stem>=8: background p_bg 0.00210, editor MH 2.969, calibrator MH 1.660
    EXCESS stem>=8 calls = (2.969-1.660) x 0.00210 x 83,530 = 230
Those 230 excess hairpin calls must come from the editor's OWN M mutations, which therefore
need a hairpin rate h8 satisfying M x (h8 - h_endo) = 230, with h_endo = 0.00349:
        M        required h8    x background
        453         0.5104         243.0x
      1,813         0.1301          62.0x
     10,000         0.0264          12.6x
     50,000         0.0081           3.8x
ENDOGENOUS A3A HAIRPIN SELECTIVITY TOPS OUT AROUND 2-3x.
AND THE BURDEN ENDPOINT CAPS M: ratio 0.92-0.98 against a 12% clone floor means M <= ~10,000.
At every M the burden endpoint ALLOWS, the required selectivity is 12-243x -- implausible for
a variant ENGINEERED to be hypoactive. At the only M where selectivity is plausible (~50,000),
M would be 60% of all calls, which the burden endpoint EXCLUDES.
THE TWO CONSTRAINTS INTERSECT IN AN EMPTY SET. The +0.318/+0.396 cannot be an editor mixture
effect at any mixing fraction this data permits. The parsimonious reading is a CLONE-LEVEL
property of the call set -- exactly what VA's sign disagreement (section 37) already said.
I did not find this argument myself. It came from the panel and I verified the arithmetic.

### SECOND ERROR, MINE, AND IT INVALIDATES EVERY SMALL p I QUOTED
Site-level permutation and bootstrap nulls treat SITES as exchangeable replicates. THE
EXCHANGEABLE UNIT IS THE CLONE. With 2 editor clones vs 2 calibrator clones there are 3
distinct clone permutations, so THE MINIMUM ATTAINABLE p IS 0.33. Every "z=5.0", "p_ed 0.0000"
and permutation z I reported for an ARM-level claim is measuring within-clone site count, not
an arm effect. Those p-values are struck.

### THIRD, AND IT IS A DESIGN-LEVEL POINT I SHOULD HAVE SEEN
The editor is expressed TRANSIENTLY, BEFORE single-cell isolation. Its mutations are therefore
CLONAL (VAF ~0.5). Our endpoint sits at VAF < 0.15, which is 91-94% of calls -- post-bottleneck
mutation plus 2-of-25-read error. WE MEASURED A STRATUM IN WHICH THE EDITOR WAS NO LONGER
ACTIVE. That is an identification failure, not a power failure, and it explains the whole
picture at once: no burden excess, symmetric call sets, and a compositional shift that behaves
like a clone property because that is what it is.

### CORRECTIONS TO MY OWN FRAMING, FROM THE PANEL, ACCEPTED
 - "83,530 vs 84,280 -- THE SAME NUMBER" OVERSTATES. Calibrator-specific counts range
   79,654-90,757 depending on comparator. The honest claim is "indistinguishable within the
   noise floor", not "identical".
 - "6x smaller than the noise" is 5.5x at 2.2% and 24x at 0.5%. Quote the range.
 - my 0.075 and 0.237 are |differences| of two clones, not SDs. E|X1-X2| = 1.128 sigma, so
   sigma = 0.066 and 0.210.

### THE REGIONAL / EDITING-INDEX PROPOSAL IS STRICTLY DOMINATED
The argument I should have made and did not: GENOME-WIDE SUMMATION IS MAXIMAL AGGREGATION, AND
THAT IS THE BURDEN ENDPOINT, WHICH IS NULL. Any partition into regions can beat the full
aggregate only via weights correlated with true effect density -- and estimating those weights
IS the prediction problem that returns 1.748x and loses to a single integer. Plus the measured
regional reliability of 0.033 at 2Mb attenuates any downstream correlation by sqrt(0.033)=0.18.
Regional aggregation is the full-aggregate estimator plus a weighting that cannot be estimated.

### POWER, COMPUTED
 burden, CV 12%:  Delta 2.2% -> 467 clones/arm;  1.0% -> 2,261;  0.5% -> 9,043
 compositional:   sigma 0.066 -> n~3;  sigma 0.210 -> n~7-8;  with variance overhead 10-12/arm
 CHEAPEST DECISIVE EXPERIMENT: n=5 clones each of Y130G and VA (~10 WGS) settles which
 variance floor is real -- an F-test on a 10x ratio, F_0.05(4,4)=6.39.
