# Data provenance audit (2026-05-23) — coordinate/version issues found

Repeated coordinate inconsistencies in the inherited site tables undermine deeper validation.
Each was found by a genome/motif sanity check (genome[pos] base + TpC spectrum).

| Source / file | Issue | Evidence | Status |
|---|---|---|---|
| **Lei_2021** (v2) | +2 bp coordinate offset | motif reads CCC not TpC; −2 shift restores TpC 0.29→0.60 | dropped for v1 |
| **v3 `dna_offtarget_sites.parquet`** (has `edit_rate`) | systematic coordinate offset | only 40% of Doman sites have C/G at pos (78% within ±1bp) | rates UNUSABLE until re-aligned |
| **Doman/McGrath** (v2) | strand encoding: minus-strand C-edits stored as `ref=C` (genome shows G) | 100% have C/G at pos, but only ~53% have C (rest G) | RECOVERABLE — orient by genome base |
| **BE4_clone1 WGS** | germline/clonal dominated | de-novo C→T 98.6% VAF=1.0, 25% TpC (no APOBEC enrichment) | not a usable continuous gold |

## Implications
- **v2 positions are correct** (100% C/G at pos) — the f-preference result (0.82 transfer) is valid;
  orienting by genome base (not stored strand) recovers ~2x more usable Doman/McGrath sites.
- **Continuous-rate validation (F2) is blocked**: the only per-site rates (v3) have offset coordinates.
  Fix = re-align v3 to v2 coords (find the offset, like Lei) OR re-derive rates on v2 coords.
- **No clean continuous DNA gold** exists in current data (v3 rates broken; WGS germline; v2 has no rates).

## Recommended canonical cleanup (gating for deeper validation)
Build ONE canonical site table: v2 coordinates + genome-based strand resolution + motif-QC per source
(drop/fix sources failing TpC) + rates re-aligned from v3 (offset-corrected) where available.
