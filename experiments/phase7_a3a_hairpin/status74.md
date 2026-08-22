
## ############################################################################
## QA TICK ~16:35 UTC — STRAND CHECK, motivated by BUG 4. I found a real asymmetry, chased
## it with the right control, and it resolves CONSERVATIVELY.
## ############################################################################

Bug 4 was s6_pileup counting only uppercase T, dropping half the reads STRAND-ASYMMETRICALLY.
That is fixed, but its existence is why the headline signal deserves a strand audit.

### THE SIGNAL IS ON BOTH STRANDS
    strand        n_elig      p_bg   editor   n_hp    calib   n_hp   ed-cal
    0 (plus)   107,555,404  0.01669   1.297    902    1.028    722   +0.269
    1 (minus)  107,455,365  0.01671   1.403    981    1.015    716   +0.387
Both clear the +0.11 bar with ample n. Asymmetry 0.119 -- the minus strand is stronger, which
is worth noting but is not the effect appearing on one strand only.

### THEN I FOUND SOMETHING THAT LOOKED BAD
alt_fwd/alt_rev inside the editor-specific sites, split by the C's genomic strand:
    editor      1.090 (C on plus) / 0.888 (C on minus)   spread 0.202
    calibrator  1.028              / 0.938               spread 0.090
A mirror-image read-orientation skew tied to the reference base is EXACTLY bug 4's family, and
the editor's is 2.2x the calibrator's. Shared in direction but not in magnitude.

### THE DECISIVE CONTROL: split it by HAIRPIN status
    arm   stratum        plus   minus   spread   increment at hairpins
    ed    non-hairpin   1.088   0.888   0.200
    ed    HAIRPIN       1.193   0.896   0.297        +0.097
    cal   non-hairpin   1.026   0.940   0.086
    cal   HAIRPIN       1.130   0.832   0.298        +0.212

 1. HAIRPIN SITES CARRY MORE SKEW IN BOTH ARMS, and the hairpin-stratum spreads are
    ESSENTIALLY IDENTICAL -- 0.297 vs 0.298. The elevated skew at hairpins is a PROPERTY OF
    HAIRPIN SITES, not something the editor does to them.
 2. The editor DOES carry more baseline skew -- 0.200 vs 0.086 in the NON-hairpin stratum,
    where hairpins play no role. That is a real sample-level property and I am recording it
    rather than explaining it away. It is not composition.
 3. BUT THE DIRECTION IS CONSERVATIVE. The increment from non-hairpin to hairpin is +0.097
    for the editor and +0.212 for the CALIBRATOR. If orientation-dependent calling inflated
    hairpin calls, it would inflate the CALIBRATOR MORE -- which works AGAINST the observed
    enrichment, not for it.

The concern resolves, and it resolves in the direction that does NOT help the headline stand.
That is the only kind of resolution worth trusting.

### STANDING POSITION UNCHANGED
Still EXPLORATORY. Still no editor claim. The inversion against A3A-Y130F is still unexplained.
Checks passed now also include strand stratification and a read-orientation audit with its
hairpin-split control. VA -- the other haA3A variant -- is 0.8h and 3.8h out and is the result
that would actually move this.
