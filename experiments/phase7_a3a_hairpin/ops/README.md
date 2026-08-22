# Operational drivers

These run the pipeline unattended. They lived only on the nodes until now, which meant a
lost instance took them with it.

| file | node | what it does |
|---|---|---|
| `progress.sh` | both | TRUE alignment progress: live fd offsets ÷ real file size |
| `auto_advance2.sh` | ai-chem | fires the editor protocol when the clean calibrator's counts complete |
| `auto_advance_B2.sh` | ai-chem2 | fires the 3-clone cross-family analysis when the Lj-BE clones complete |
| `keeper.sh` | ai-chem | hardlinks D10A BAMs so a stale `KEEP_PATTERN` in the driver cannot delete the calibrator |
| `s6_driver_v2.sh` | ai-chem | pileup queue |

## Reading progress.sh output

    SAMPLE   MATE1   MATE2   CUR_READS   IMPLIED   ETA_AVG   ETA_NOW

- **MATE1/MATE2** are authoritative. Both mates should agree within ~0.1%; a divergence means
  one reader is stuck.
- **CUR_READS double-counts any sample re-aligned with the same thread count.** These logs are
  append-only and two generations with the same `-K` are indistinguishable by batch size. No
  batch-size rule can separate them.
- **IMPLIED** = CUR_READS extrapolated by the file fraction. It exists to expose exactly that:
  compare it across samples with similar file sizes. D10A-clone6 reads 2618M against a true
  ~991M for its 9.3 GB file, because it was re-aligned after a preemption.
- **ETA_AVG** is a lifetime average and cannot see a slowdown. **ETA_NOW** is computed between
  the last two invocations via `/tmp/progress_state.txt`; it prints `-` on the first run after
  a restart rather than guessing.

## Two rules these encode

**Never edit a running bash script.** bash re-reads by byte offset; an edit mid-run either
does nothing or executes garbage. Write a new file. This cost one bad run.

**Never `pkill -f <pattern>` where the pattern appears in your own command line.** It matches
the shell you are typing in. This killed two SSH sessions. Use the bracket form `"[c]hain.sh"`
or kill by PID.
