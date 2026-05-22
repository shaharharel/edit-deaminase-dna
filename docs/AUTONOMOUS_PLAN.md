# Autonomous execution plan (post-rethink, 2026-05-23)

Context: the deaminase-preference (f) model is held-out-validated cross-cell-type (AUROC 0.82,
sequence-only ±10bp probe). Chromatin (g) is unvalidated; WGS continuous gold failed (clone germline).
Assets in hand: train_table_v1 (handcrafted f + ENCODE g, 75,924 rows); HyenaDNA union embeddings
(256d, 120,160 sites, /tmp/hyenadna_union.pt); seq_1kb.parquet; per-site rates (Doman/Richter/Yu);
Marquart BE-DICT efficiency libraries (MOESM3 ABE / MOESM4 CBE).

## Validation modality map (why each matters)
- **BE-DICT / efficiency libraries** → precise, high-N, per-editor gold for the **f (deaminase chemistry)**.
- **Per-site off-target rates** → f × local-exposure at real off-target sites (medium power).
- **WGS (bulk deep)** → only modality for **g (genome-wide location/burden)**; acquisition-gated; clonal fails.

## F1 — Full-model vs simple-probe comparison  [AUTONOMOUS, data in hand]
Compare feature stacks on the SAME held-out transfer test (train Doman → test McGrath, rAPOBEC1):
- (a) seq ±10bp one-hot (the 0.82 probe), (b) handcrafted-81, (c) HyenaDNA-256, (d) all combined.
- Report AUROC per feature set. NOTE: union-pool negatives are trinuc-matched (Strategy B) → removes the
  motif → "beyond-trinuc" operating point. For an apples-to-apples vs 0.82 (random negatives), also:
- F1b: re-extract handcrafted + re-embed HyenaDNA for a RANDOM genomic negative set, rerun on random negs
  (the true 0.82-comparable comparison). HyenaDNA re-embed via scripts/embed/compute_hyenadna.py.
- Deliverable: feature-ablation table; does richer representation beat the ±10bp probe.

## F2 — Continuous-rate validation = the f gold  [AUTONOMOUS, data in hand]
- Parse Marquart BE-DICT MOESM3 (ABE) / MOESM4 (CBE): context → editing rate per editor.
- Test: does the preference model predict BE-DICT editing rates (Spearman/R²), per editor? Disentangle
  window-position from sequence motif (control for position-in-window).
- Test: predicted preference vs measured per-site VAF (Doman/Richter/Yu), guide-independent subset.
- Deliverable: continuous-validation table (the well-powered f claim).

## F3 — Honest deliverables  [AUTONOMOUS]
- Per-editor deaminase-preference ranking (from BE-DICT + our model).
- Intrinsic-susceptibility per-gene prior + COSMIC/ClinGen Tier A/B/C/D stratification (download lists).
- Update MODEL_CARD with F1/F2 results + calibrate where a continuous gold exists.

## G — Bulk-WGS gold acquisition  [GATED, research + download]
- Identify bulk treated+control deep-WGS BE datasets (Buchumenski 2021 supp datasets; recent studies).
- Multi-clone/multi-control filtering recipe (à la Doman) to extract clean APOBEC edits if only clones exist.
- This unlocks the g/regional/absolute-burden validation; not blocking F1–F3.

## Guardrails (carry into every step; from QA)
- Hold negatives fixed when comparing feature sets; report which negative regime.
- TpC-spectrum QC gate on any WGS-derived signal before trusting it.
- Chromosome-block holdout + permutation null for any correlation; report block-bootstrap CIs.
- Never match negatives on model-own features; honest resolution (calibrate only ≥1 Mb).
- Run ml-code-reviewer + scientific-analyst QA on each new stream.
