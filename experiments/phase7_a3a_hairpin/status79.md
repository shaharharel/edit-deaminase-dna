
## ############################################################################
## 2026-08-22 ~18:00 UTC — "IS IT LEARNABLE AT ALL?" ANSWERED. Yes, at 4.422x --
## and the deaminase-free control is learnable at 2.350x, which is the whole story.
## ############################################################################

Trained ON THE EDITOR'S OWN SITES for the first time. Held-out CHROMOSOME, negatives matched
on trinucleotide AND strand, base rate 0.0909 so the ceiling is 11.00x, random baseline beside
every number. This is the CEILING of what this data supports -- no transfer, no borrowing.

    topK    ED enr  ED rand   CAL enr  CAL rand   ED/CAL   sigma
    0.1%     4.422    0.935     2.350     1.080     1.88     7.3
    1.0%     1.898    0.928     1.538     1.026     1.23     5.6
    5.0%     1.359    0.995     1.238     1.002     1.10     4.9

### 1. YES, IT IS LEARNABLE -- 4.422x at top 0.1% against a random 0.935
That is a real, well-powered signal and it is the first time this project has trained on the
editor data itself rather than trying to transfer into it.

### 2. BUT THE DEAMINASE-FREE CALIBRATOR IS ALSO LEARNABLE, AT 2.350x
nCas9 has no deaminase. Whatever the model finds in its sites is ENDOGENOUS somatic mutation
structure. So MOST OF WHAT IS LEARNABLE HERE IS NOT THE EDITOR. This is the same lesson the
whole session has taught in five different ways, now measured directly.

### 3. THE EDITOR INCREMENT IS REAL BUT MODEST: 1.88x over the endogenous baseline, 7.3 sigma
And it decays with panel size exactly like every other signal in this project:
1.88x -> 1.23x -> 1.10x across 0.1% -> 1% -> 5%.

### 4. AND IN-DOMAIN TRAINING DOES NOT FIX THE CAPTURE PROBLEM
    top 0.1%:    918 sites hold   369 of 83,530 editor mutations =  0.44% captured
    top 1.0%:  9,188 sites hold 1,585 of 83,530                  =  1.90% captured
    top 5.0%: 45,941 sites hold 5,675 of 83,530                  =  6.79% captured
At the depth where the model is strongest (4.42x) it holds 0.44% of the editor's mutations.
ENRICHMENT AND CAPTURE ARE IN DIFFERENT PLACES NO MATTER WHAT YOU TRAIN ON. That was true for
the PCAWG model and it is equally true for a model trained on the editor's own data, which
means it is a property of the DATA, not of the training source.

### WHAT THIS SETTLES
The question was never "can a cancer model transfer" -- section 32 already showed `stem` alone
beats the cancer model. The question was whether the editor's off-targets are structurally
predictable at all. THEY ARE, at 4.42x, of which 2.35x is endogenous and 1.88x is the editor.
The obstacle to a usable safety panel is not the model, the architecture, or the training set.
IT IS THAT THE EDITOR'S MUTATIONS ARE DIFFUSE: the top 0.1% of the most predictable sites in
the genome contains one in every 226 of them.

### HEALTH
ai-chem TERMINATED by my own stop, queue complete. ai-chem2 3 bwa, 0 failed, 199G,
VA-clone2 77.7% (1.4h), YE1-clone1 10.2% (5.9h).
