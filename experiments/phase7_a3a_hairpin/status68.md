
## ############################################################################
## QA TICK ~13:15 UTC — I APPLIED THE clone7 TEST TO THE haA3A POSITIVE. IT SURVIVES.
## (node A is STOPPED, so this is recorded in the repo and on node B, not /mnt/data/a3a)
## ############################################################################

Half an hour ago I disqualified A3A-Y130F-clone7 for carrying 62x its siblings' specific
sites, and threw out the 3-way result it produced. The obvious next question is whether the
Y130G positive in section 25 has the same disease. I had gated the CALIBRATOR but never the
EDITOR CLONES that produced the effect.

### THE Y130G CLONES BOTH QUALIFY -- on the same gate that rejected clone7
    Y130G-clone1   790,444 alt>=1   46.5% sub-VAF-0.05   cov 21.77   0.92x peer median  EXIT 0
    Y130G-clone2   859,960 alt>=1   55.8% sub-VAF-0.05   cov 24.52   1.01x peer median  EXIT 0
Compare what failed: clone7 14,818,183 / 91.3% / 7.72x, background 16,960,881 / 90.3% / 9.63x.
The Y130G clones are in the clean band with nCas9 and Parent, not near the rejects.

### THE SUSPICIOUS COINCIDENCE, CHECKED RATHER THAN ASSUMED
The haA3A output paired up too neatly to trust at face value:
    Y130G-clone2 90,513 vs nCas9-clone1 90,517   -- 4 apart
    Y130G-clone1 80,028 vs nCas9-clone2 80,053   -- 25 apart
Two independent pairs matching to under 0.03% is usually a bug -- the same number printed
twice, or crossed labels. I RECOMPUTED the specific-site counts myself from the counts files
rather than trusting the script:
    pair                              jointly eligible   ed specific   cal specific   ratio
    Y130G-clone2 vs nCas9-clone1          215,010,769         83,530        84,280   0.991
    Y130G-clone1 vs nCas9-clone2          214,268,113         74,077        73,247   1.011
    Y130G-clone1 vs nCas9-clone1          214,044,700         74,253        84,517   0.879
NOT A BUG. The third row proves the script is not printing one number twice -- that pairing
comes out clearly unequal at 0.879. The near-equality in the first two rows is a REAL property
of these samples: they carry near-identical specific-site burdens.
(My absolute counts differ from the script's because my eligibility definition is not
byte-identical to s7c's; the RATIOS are the point and they confirm the pairing.)

### WHAT THAT MEANS FOR SECTION 25
The haA3A comparison is BURDEN-MATCHED TO WITHIN 1% -- the cleanest editor-vs-control burden
matching anywhere in this project. A3A-Y130F was 8.7x off. Every check I can construct now
passes on this arm: gate A0 on both editor clones, burden matched 0.991/1.011, GC-shift
asymmetry 0.4pp, measured clone-luck floor +0.075 against observed +0.318/+0.396, true-null
validation +0.037/+0.003/-0.057, and monotone scaling with stem length.

### AND I STILL HOLD THE FRAMING
The pre-registration declared this EXPLORATORY and expected no enrichment. It remains
exploratory. I have not converted it into a confirmatory finding and I am not claiming haA3A
has hairpin off-target activity. What I can say is narrower and it is what the data support:
the positive has now passed every control that killed other results this session, and THE
BIOLOGICAL INVERSION AGAINST A3A-Y130F IS STILL UNEXPLAINED. Until it is, no editor claim.
