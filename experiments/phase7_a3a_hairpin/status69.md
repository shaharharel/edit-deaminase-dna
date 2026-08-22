
## ############################################################################
## 2026-08-22 ~13:30 UTC — A CONCRETE, QUANTIFIED LEAD ON THE INVERSION. The two arms'
## analysis sets are NOT the same population, and the divergence is in one filter step.
## ############################################################################

Sections 25 / 25.1 said no editor claim can be made until the biological inversion is
explained -- A3A-Y130F (active) depleted, Y130G (engineered-safe) enriched. This is the first
concrete lead, computed from numbers already on record since node A is stopped.

### WHERE THE A3A-Y130F SITES GO
    sample                  cov     alt>=1    alt>=2  specific  spec/alt>=2
    A3A-Y130F-clone2      25.52    864,617   184,003    19,841     0.108
    A3A-Y130F-clone5      26.22    861,941   185,005    19,694     0.106
    nCas9-clone1          24.25    854,914   256,302    90,757     0.354
    nCas9-clone2          24.25    780,594   245,984    80,053     0.325
    D10A-clone6           39.97  5,837,065   346,239   171,632     0.496

AT alt>=2 the A3A-Y130F clones hold 72% of what nCas9 holds (184k vs 256k) -- a modest gap
that coverage and clone burden could account for.
AFTER the specificity filter they hold 22% (19.8k vs 90.8k). THE GAP TRIPLES IN THAT ONE
STEP, and it is precisely the step that defines the analysis set every MH OR is computed on.
    A3A-Y130F retains 10.6-10.8% of its alt>=2 sites
    nCas9      retains 32.5-35.4%
    D10A-clone6 retains 49.6%

### AND IT IS NOT THE CALIBRATOR VETO, WHICH IS THE OBVIOUS SUSPECT
The specificity filter drops a site if the CALIBRATOR calls it. But the calibrator's own
alt>=2 rate is ~256,302 / 215M eligible = 0.12%. A 0.12% veto cannot mechanically remove 89%
of anything. The loss must sit in the Parent germline mask or in joint eligibility, and
whatever it is, IT HITS THE A3A-Y130F CLONES ABOUT 2.5x HARDER THAN IT HITS nCas9.

### WHY THIS MATTERS MORE THAN IT LOOKS
Every editor number in sections 21, 21.2, 23.1 and 25 is computed on the "specific" set. If
that set is assembled by a filter that removes 89% of one arm and 65% of the other, the two
arms are not measuring the same population of sites, and a difference between their hairpin
fractions need not be an editor difference at all. This is a candidate explanation for the
inversion that does not require either editor to behave strangely.
IT IS A LEAD, NOT A FINDING. I have not shown the filter is responsible, only that it is where
the arms diverge and that the obvious mechanism (calibrator veto) cannot account for it.

### NEXT ACTION, AND IT NEEDS NODE A
Per-filter-stage counts -- how many sites survive joint eligibility, then the Parent mask, then
the calibrator veto, separately for each arm. Those need the counts files on node A, which is
stopped. THIS IS THE FIRST THING TO RUN IF NODE A IS RESTARTED, ahead of any new arm.

### HEALTH
ai-chem  TERMINATED -- stopped deliberately last tick after its queue completed and a clean
         co-tenant check, NOT a preemption. No restart action taken or needed.
ai-chem2 3 bwa, 0 failed, 207G, 87.7% us. VA-clone1 54.4% (2.5h), VA-clone2 12.1% (5.5h),
         YE1 x2 queued. VA is the replication that matters most: it is the OTHER haA3A
         variant, so it independently tests the section 25 positive.
