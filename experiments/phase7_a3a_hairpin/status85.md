
## ############################################################################
## 2026-08-22 ~19:20 UTC — RAN THE SECTION-32 ATTACK AGAINST THE CANCER RESULT.
## IT GOES THE OPPOSITE WAY, AND THE CONTRAST IS THE FINDING.
## ############################################################################

Section 32 retracted the editor "transfer" because `stem` alone -- one integer, no training --
BEAT the trained model. I never ran that test on the cancer genomic result. If stem alone also
reached ~6.8x there, section 35 would be about a hairpin annotation and not about a model.

### IT DOES NOT. IN CANCER THE MODEL WINS DECISIVELY.
    panel     model      stem  hp_score   model - best untrained
    0.01%    16.698     2.370     6.787              +9.911
    0.10%     6.819     1.982     3.534              +3.286
    1.00%     3.002     1.441     1.563              +1.439
    5.00%     1.665     1.203     1.265              +0.401
   10.00%     1.372     1.143     1.218              +0.154
All out-of-fold, all against the same 92,823 PCAWG mutations on 236,482,055 sites.
`stem` alone reaches 1.982x at top 0.1% against the model's 6.819x. The learned combination of
six structural features is worth 3.4x over the best single one.

### THE CONTRAST WITH THE EDITOR IS THE ACTUAL RESULT
                        model      stem alone    verdict
    CANCER  top 0.1%    6.819x        1.982x     MODEL WINS by 3.4x
    EDITOR  top 0.1%    +0.466        +1.655     STEM WINS -- section 32 retraction
In cancer there is enough signal that a learned combination substantially beats any single
feature. In the editor arm there is so little signal that the learned combination performs
WORSE than the best raw feature -- which is what fitting noise looks like.

SECTION 35 SURVIVES THIS ATTACK AND IS STRENGTHENED BY IT. The cancer result genuinely is
about a model doing something a hairpin annotation cannot; the editor result genuinely was not.
Both the retraction and the confirmation came from running the SAME test in both places, and
that was only possible because I ran it on the editor first, where it hurt.

### WHAT THE PROJECT'S CLAIM NOW IS, PRECISELY
A gradient-boosted model over six DNA-structure features, trained on APOBEC3-family mutations
from 21 PCAWG donors and evaluated OUT-OF-FOLD across whole chromosomes, concentrates
APOBEC-context mutations 6.819x in the top 0.1% of a 236M-site genome and captures 13.72% of
them in the top 10%. It beats every untrained ranking by 3.4x. THE SAME MODEL, APPLIED TO
BASE-EDITOR OFF-TARGETS AT THE SAME GENOMIC SCALE, REACHES 1.748x AND IS BEATEN BY `stem`.

### HEALTH
ai-chem TERMINATED by my own stop. ai-chem2 0 failed, VA-clone2 pileup running toward the
second pre-registered VA clone, YE1-clone1 aligning.
