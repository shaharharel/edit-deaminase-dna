
## ############################################################################
## QA TICK ~19:00 UTC — SWEPT THE RECORD FOR THE DENOMINATOR DEFECT. The document is
## clean; THE PUBLISHED ARTIFACT WAS NOT. Fixed the deliverable.
## ############################################################################

### THE SWEEP OF HANDOFF -- CLEAN, but my scanner tested the wrong thing
167 enrichment-like values; 5 lacked a nearby base-rate / denominator / baseline reference,
and all 5 are cross-references inside the section-20 disambiguation paragraph -- false
positives. So 162/167 carry their context.
BUT THAT IS NOT THE DEFECT THAT OCCURRED. The section-33 error was not a MISSING base rate --
0.0909 was stated. It was two STATED base rates being compared across sections as if
commensurable. My scanner cannot see that, and I am recording that limitation rather than
claiming the sweep proves more than it does.

### THE REAL TEST FOUND IT IN THE PUBLISHED ARTIFACT
    HEK293T panel rows   base rate 83,530 / 215,010,769 = 0.000388   1 in 2,574
    PCAWG in-sample      base rate 83,999 /     923,989 = 0.0909     1 in 11
    THE TWO BASE RATES DIFFER BY 234x.
The artifact's standfirst read "predicts held-out tumour chromosomes at 5.369x ... ranks the
editor's off-target sites at 0.95-1.00x", and the table carried 5.369x as a "positive control"
directly beneath the HEK293T rows. Both figures are legitimate enrichments over their OWN base
rate, so "clear skill vs none" is right -- but a reader differences them and gets a 5.6x drop,
and part of that gap is nothing but the denominator.
THIS IS THE ONE PLACE IT MUST NOT BE AMBIGUOUS, because the artifact is the thing that leaves
this machine.

### FIXED IN THE DELIVERABLE
 - standfirst now states both base rates, the 234x ratio, and instructs the reader to take the
   numbers as "clear skill versus none", NOT as a 5.6x drop.
 - table caption now names the HEK293T base rate (1 in 2,574) at the rows, and labels the
   in-sample control as being on the TRAINING SET's 1-in-11 denominator, "included to show the
   model has skill on its own data, not to be differenced against the rows above".
Republished to the same URL. Tag balance verified before publishing.

### WHAT I AM NOT CLAIMING
That the sweep proves the rest of the record is free of this defect. It proves every
enrichment states its own base rate. It does NOT prove that no two of them are compared
across scales somewhere I have not looked. The only defence against that is the rule added
last tick -- state the denominator with the value, every time -- and it is a habit, not a test.

### INTEGRITY -- node B clean, node A stopped with its queue complete
