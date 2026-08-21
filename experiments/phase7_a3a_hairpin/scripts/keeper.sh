#!/bin/bash
# Preserve the CALIBRATOR BAMs that the running s6 driver would delete.
#
# The driver was started with KEEP_PATTERN='A3A-Y130F', so P66-D10A-* BAMs are freed
# after their 23/23 counts. D10A is the deaminase-free CALIBRATOR for the editor test -
# half of the comparison - and regenerating one is ~13 h. Retaining the editor arm but
# not its calibrator is asymmetric: if the editor shows elevated burden, the first
# forensic question is whether the CALIBRATOR has a mapping artefact, and that needs
# its alignments.
#
# Restarting the driver to pick up a corrected pattern would interrupt an in-flight
# pileup, so instead: hardlink every D10A bam on sight. A hardlink costs zero bytes
# (same inode); the driver's `rm -f "$bam"` removes only its own name and the data
# survives under KEEP_.
cd /mnt/data/a3a/wgs || exit 1
echo "[$(date '+%F %T')] keeper started; hardlinking P66-D10A-* bam/bai on sight"
while true; do
  for f in P66-D10A-*.bam P66-D10A-*.bam.bai; do
    [ -e "$f" ] || continue
    case "$f" in KEEP_*) continue ;; esac
    if [ ! -e "KEEP_$f" ]; then
      ln -f "$f" "KEEP_$f" && echo "[$(date '+%F %T')] hardlinked KEEP_$f"
    fi
  done
  sleep 120
done
