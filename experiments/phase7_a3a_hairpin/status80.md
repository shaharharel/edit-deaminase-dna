
## ############################################################################
## CORRECTION ~18:30 UTC — SECTION 33 QUOTED TWO ENRICHMENTS ON DIFFERENT DENOMINATORS
## WITHOUT SAYING SO, AND THE CAPTURE FIGURE I GAVE WAS 2.6x TOO FLATTERING.
## Found because the user asked what n was behind the 4.422x.
## ############################################################################

### THE TWO NUMBERS ARE NOT ON THE SAME SCALE
    LEARNABILITY TEST (s13)
      denominator  the CONSTRUCTED trainset: 918,830 rows = 83,530 pos + 10x matched neg
      top 0.1%     918 sites, ~369 editor positives   <- n is NOT low, it is well powered
      base rate    0.0909 (1 in 11) BY CONSTRUCTION, because I downsampled negatives 10:1

    GENOMIC PANEL TEST (s12)
      denominator  ALL jointly-eligible sites: 215,010,769
      base rate    0.000388 = 1 in 2,574   <- the REAL genomic rate
      top 0.1%     215,010 sites, 146 mutations captured, 1.748x

4.422x is enrichment over a base rate I CREATED. It answers "can a model tell an editor site
from a MATCHED non-site" -- yes, and n=369 is ample. 1.748x is enrichment over the real
genomic rate and answers "can a model FIND editor sites in the genome" -- barely.
The trainset negatives are a 257x downsample of the genome, so the learnability number
CANNOT be read as genomic performance. Section 33 put both in the same discussion without
flagging that, which invites exactly the wrong reading.

### AND THE CAPTURE FIGURE WAS WRONG IN THE FLATTERING DIRECTION
Section 33 said "0.44% captured". That was 369/83,530 -- computed from the TRAINSET panel.
The GENOMIC panel at the same 0.1% depth captures 146/83,530 = 0.17%.
    genomic top 0.1%    215,010 sites hold   146 of 83,530 = 0.17%
    genomic top 1.0%  2,150,107 sites hold   957 of 83,530 = 1.15%
    genomic top 5.0% 10,750,538 sites hold 4,114 of 83,530 = 4.93%
THE CORRECT NUMBER IS 2.6x WORSE THAN THE ONE I REPORTED. The conclusion is unchanged --
both are far too small for a panel -- but I quoted the better of two numbers that came from
different denominators, and I did not notice until asked what n was.

### WHY THIS MATTERS BEYOND THE ARITHMETIC
This is the same bug family the whole session has been about, in my own reporting rather than
in code: a quantity correct in its own context (enrichment over a 1:11 constructed base) used
downstream as if it were universal (enrichment over the genome). It is the seventh instance
and the first one that reached a number I gave the user.
RULE ADDED: every enrichment quoted from now on states its DENOMINATOR and its BASE RATE, not
just its value. "4.422x" is meaningless without "over a 0.0909 base rate on a constructed set".
