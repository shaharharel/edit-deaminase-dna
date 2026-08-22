
## ############################################################################
## QA TICK ~17:35 UTC — I ATTACKED SECTION 31'S FRAMING AND IT DOES NOT SURVIVE.
## RETRACTING "the cancer-trained model transfers". `stem` alone does it better.
## ############################################################################

### THE ATTACK
Y130G's DEFINING property is hairpin excess (+0.318/+0.396). The PCAWG model ranks by hairpin
geometry. So a hairpin-ranked panel will of course contain more Y130G sites. The "transfer"
may be nothing but the hairpin enrichment re-expressed as a ranking, with the cancer training
contributing nothing. Test: race the trained model against rankings that need NO training and
NO cancer data.

### THE RESULT, AND IT IS DECISIVE
    ranking          panel   n_ed    ed x  n_cal   cal x      gap
    PCAWG model      0.01%     35   4.190     12   1.424   +2.766
    stem ALONE       0.01%     31   3.711      8   0.949   +2.762
    hp_score ALONE   0.01%     30   3.592     14   1.661   +1.930
    PCAWG model      0.10%    146   1.748    108   1.281   +0.466
    stem ALONE       0.10%    275   3.292    138   1.637   +1.655
    hp_score ALONE   0.10%    225   2.694    131   1.554   +1.139
    PCAWG model      1.00%    956   1.144    927   1.100   +0.045
    stem ALONE       1.00%   1331   1.593    903   1.071   +0.522

    PCAWG minus best-untrained:  +0.004 at 0.01%,  -1.188 at 0.10%,  -0.477 at 1.00%

`stem` IS ONE INTEGER READ STRAIGHT OUT OF THE UNIVERSE FILE. No training. No cancer data.
It TIES the PCAWG model at 0.01% and BEATS IT DECISIVELY at 0.1% and 1% -- the panel sizes
that matter. THE CANCER TRAINING CONTRIBUTES NOTHING AND ACTIVELY DEGRADES THE RANKING.

### WHAT I AM RETRACTING
Section 31's framing -- "a cancer-trained hairpin model ranks base-editor off-target sites" --
is UNSUPPORTED. The correct statement is the one already known from the endpoint: hairpin-rich
sites are enriched among Y130G's calls. The panel was that same fact wearing a ranking.
WHAT SURVIVES: the CONDITIONAL pattern is still real and still interesting -- ranking by
hairpin structure separates editor from calibrator in Y130G (+2.76) and not in VA (-0.15) or
A3A-Y130F (-0.48). But it is a statement about HAIRPIN GEOMETRY, not about cancer training,
and the best ranker is the simplest possible feature.
NOTE FOR THE RECORD: if a usable panel is ever built here, it should be built from `stem`,
not from a trained model. stem-alone at 0.1% gives 3.292x with n=275 -- better than anything
the model produced at any depth.

### AND THE NUMBERS BEHIND THE EDITOR DATA, ASKED FOR DIRECTLY
WGS, clonal HEK293T, hg19, ~236M TCW universe sites per sample.
    PRJNA1042830  Y130G x2, VA x2, YE1 x2, nCas9 x2, Parent    cov 21.8-27.9
    PRJNA1006866  A3A-Y130F x3, D10A x3, background            cov 25.0-40.0
Editor-specific positives: Y130G-c2 83,530 / Y130G-c1 80,028 / VA-c1 94,361 /
A3A-Y130F-c2 19,841 / A3A-Y130F-c5 19,694.
HOW SURE ARE WE THEY ARE EDITOR? AT THE SITE LEVEL, NOT AT ALL:
    Y130G-clone2 editor-specific     83,530
    nCas9-clone1 calibrator-specific 84,280   <- DEAMINASE-FREE, so entirely endogenous
The deaminase-free clone yields essentially the SAME NUMBER. The burden endpoint agrees --
0.922/0.944/0.959/0.978, at the clone floor, NO excess. There is no count excess to attribute
to the editor anywhere in this dataset. The only editor-attributable signal is COMPOSITIONAL
(the hairpin fraction), never the number of mutations.

### LAUNCHED: IS IT LEARNABLE AT ALL?
s13_learnable.py trains ON THE EDITOR'S OWN SITES -- held-out chromosome, trinuc+strand-matched
negatives, tail enrichment with a random baseline -- and runs the IDENTICAL pipeline on the
deaminase-free calibrator's own sites. If both are equally learnable, the model is learning
somatic mutation structure, not editor activity.
