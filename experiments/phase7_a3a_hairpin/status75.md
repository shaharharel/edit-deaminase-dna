
## ############################################################################
## QA TICK ~15:35 UTC — PREVENTIVE QA. VA-clone1 lands in ~0.3h and NOTHING would have run
## it. Armed a gated driver AND pre-registered a discriminating prediction before the data.
## ############################################################################

### THE GAP, FOUND BEFORE IT COST ANYTHING
The haA3A driver was written for Y130G only and has already exited. VA-clone1 reaches 23/23
in about 18 minutes and no driver was watching for it. The clone7 driver taught this lesson
expensively: armed without gate A0, it fired on a clone carrying 62x its siblings' specific
sites and produced a 3-way result that had to be discarded.

### WHAT I ARMED
auto_advance_va.sh, verified running by the log it writes (pid 344784), processing each VA
clone as it completes. Per clone it:
  1. RUNS GATE A0 ON THE EDITOR CLONE ITSELF -- the clone7 lesson was that I gated the
     calibrator and never the editor. On failure it writes BLOCKED.md and does NOT run.
  2. Runs gate A0 on nCas9-clone1 as well.
  3. Runs s7c_editor, then GC-stratified, then COMPLEXITY-stratified, then the strand +
     read-orientation audit -- i.e. every check the Y130G positive had to pass, applied
     from the start rather than retrofitted over five ticks.

### THE PRE-REGISTRATION, WRITTEN INTO THE DRIVER BEFORE THE DATA EXISTS
The STANDING expectation is UNCHANGED: haA3A was engineered for near-background off-target,
so the pre-declared expectation is STILL NO ENRICHMENT. Y130G contradicting it does NOT
license me to quietly flip the expectation for VA, and I have not.
What I added is a DISCRIMINATING PREDICTION whose entire value is that it is recorded now:
    VA POSITIVE  (ed-cal > +0.11, near Y130G's +0.318/+0.396)
        -> the effect tracks the haA3A CLASS: two independent variants, four clones.
    VA NULL      (inside +/-0.11)
        -> Y130G's positive is Y130G-specific; the class expectation survives and Y130G
           becomes a single-variant anomaly needing its own explanation.
    VA NEGATIVE  (< -0.11, like A3A-Y130F's -0.18/-0.21)
        -> the endpoint produces variant-specific SIGNS with no biological ordering, which is
           evidence the ENDPOINT is the problem -- exactly what the pre-registration warned.
I commit to reporting whichever lands and to NOT reinterpreting the categories afterwards.

### WHY THIS IS THE RIGHT USE OF A QA TICK
Every check I ran on Y130G was retrofitted AFTER seeing a surprising number. Retrofitted
checks are worth less than pre-committed ones, however carefully done. VA is the chance to
have the checks and the prediction in place first, and it costs nothing to do it now.

### INTEGRITY -- node B clean; node A stopped deliberately with its queue complete
