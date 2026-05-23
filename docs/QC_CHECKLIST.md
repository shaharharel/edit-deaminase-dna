# Mandatory QC gates for every new BE off-target dataset (QA tracker, G0-G9)
A source feeds the g-LANDSCAPE only if it passes G1,G2,G3(with control),G4=bulk + eval G5,G7,G8,G9.
Clonal/uncontrolled sources may feed f(motif) ONLY, with caveat.

- G0 Provenance + single pinned hg38 (md5); record editor/celltype/assay/build.
- G1 Coordinate-offset gate: >=99% C/G(or A/T) at center d=0 AND d=0 = argmax of motif-spectrum scan
     (catches Lei +2 / wrong build). MANDATORY, not diagnostic.
- G2 Motif-spectrum gate: TpC(CBE)/TpA(ABE) >40% oriented by GENOME base; NO family exemption (engineered
     editors must still show coherent non-uniform spectrum).
- G3 Germline/clonal contamination: (a) treated-minus-control; (b) VAF not spiked at 0.5/1.0; (c) TpC
     enrichment vs trinuc background. BE4_clone1 failed (98.6% VAF=1.0, 25% TpC). -> scripts/data/contamination_qc.py
- G4 Bulk-vs-clonal classifier: Gini>=0.7 + r_DNase~0.4-0.5 = bulk (trust); Gini<=0.3 + r_DNase~0 = clonal (reject for landscape).
- G5 STRICT baseline {logC,logTpC,mappability,GC,gene_density}: any track claim = delta over STRICT (CI excl 0). Bake GC+gene_density into bins.
- G6 Mappability/blacklist filter before binning.
- G7 LOCO leakage: standardize on train only; per-chrom held-out.
- G8 Leave-one-SOURCE-out across ALL sources (not hand-coded pair).
- G9 Circular-shift null WITHIN chromosome; real >> null 95th.
