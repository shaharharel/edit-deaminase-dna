# Coordinate conventions — READ BEFORE JOINING ANY TWO ARRAYS

Two bugs tonight came from strand info being destroyed by pyrimidine-orientation
and then used as if it survived. A third of the same family was caught latent:
two arrays in this directory used the field name  with OPPOSITE conventions.

| file | field | convention |
|---|---|---|
| pcawg/snvs.npz, feat/a3a_trainset*.npz | pos  | **1-based** (MAF Start_position) |
| feat/universe_chr*.npz                 | pos  | **0-based** array index |
| feat/universe_chr*.npz                 | pos1 | **1-based** (added; use THIS for joins) |
| feat/universe_chr*.bed                 | -    | 0-based half-open (BED standard) — correct |

samtools mpileup / faidx emit **1-based** positions. Join against , never .

MANDATORY assertion for any script that joins these: verify the reference base at
the joined coordinate is C (strand==0) or G (strand==1) for a sample of sites,
and abort if not. An off-by-one here does not crash — it silently converts TCW
into CWN and produces a plausible but wrong answer.

Verified 2026-08-10: universe pos is 0-based (8/8 samtools faidx checks correct
as 0-based, 0/8 as 1-based).
