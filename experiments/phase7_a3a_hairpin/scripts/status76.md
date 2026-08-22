
## ############################################################################
## 2026-08-22 ~16:40 UTC — VA LANDED AND IT IS NULL. It falls in the category I
## pre-registered at 15:35, and I am reporting it as written rather than reinterpreting it.
## ############################################################################

### VA-clone1 vs nCas9-clone1, both QUALIFIED on the corrected gate (0.65x and 1.06x)
    banded MH   stem6  editor 0.985  calibrator 1.015  ->  -0.030
                stem7  editor 1.186  calibrator 1.344  ->  -0.158
                stem8  editor 1.344  calibrator 1.612  ->  -0.268
    GC-adjusted MH     editor 0.975  calibrator 1.041  ->  -0.066
    complexity  MH     editor 0.988  calibrator 1.054  ->  -0.066

### AGAINST THE PRE-REGISTRATION WRITTEN BEFORE THE DATA EXISTED
    positive (>+0.11)   -> haA3A CLASS effect
    null (inside 0.11)  -> Y130G-SPECIFIC; the class expectation survives
    negative (<-0.11)   -> the endpoint produces variant-specific signs and is itself broken
OBSERVED -0.066 ON THE GC-ADJUSTED ENDPOINT = THE NULL CATEGORY, cleanly.
THE STANDING PRE-REGISTRATION IS UPHELD FOR VA: haA3A was engineered for near-background
off-target and VA shows none. Y130G is now a SINGLE-VARIANT ANOMALY needing its own
explanation, exactly as the null branch said it would be.
I committed to not reinterpreting the categories after seeing the number. I have not.

### THE PANEL RESULT REPLICATES, WITH A PROPER PERMUTATION NULL
    panel     clone1 gap   clone2 gap   z (clone1, 200 permutations)
    0.01%       +2.620       +2.766       5.0
    0.10%       +0.459       +0.466       2.7
    1.00%       +0.098       +0.046       1.8
    5.00%       +0.060       +0.024       2.5
Null centred at ~0.000 with sd 0.525 / 0.173 / 0.050 / 0.023. The two clones agree to 0.007
at 0.1%. This fixes the two defects in my first run -- no baseline at 0.01%, and a single
shuffle instead of a null distribution.

### THE COHERENT PICTURE, AND IT IS NARROWER THAN THE PANEL NUMBER ALONE SUGGESTS
    Y130G      hairpin excess +0.318/+0.396   AND panel transfer z=5.0, replicated
    A3A-Y130F  NO hairpin excess -0.18/-0.21  AND panel ANTI-transfer -0.483
    VA         NO hairpin excess -0.066       -> panel should be null; test running
THE MODEL TRANSFERS EXACTLY WHERE HAIRPIN SIGNAL EXISTS. That is mechanistically coherent and
it is the best explanation this project has produced for why the transfer failed the first
time. It also CONFINES THE PHENOMENON TO ONE VARIANT OUT OF THREE.

### WHAT THIS DOES TO THE HEADLINE
It does NOT become "we can predict base-editor burden". The panel that shows z=5.0 captures
35 of 83,530 editor mutations -- it misses 99.96% -- and the skill is gone by the panel size
where capture becomes meaningful. Enrichment and capture are in different places.
The defensible statement is narrower and it is about MECHANISM: a cancer-trained hairpin model
ranks off-target sites in the one editor variant that shares its sequence preference, and
fails or anti-transfers in two variants that do not.
