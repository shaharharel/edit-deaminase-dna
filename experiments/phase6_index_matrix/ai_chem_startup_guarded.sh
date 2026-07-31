#!/bin/bash
# ai-chem startup: McGrath AncBE4max WGS reprocess (hg38, /mnt/data, bio env). Resume on preemption/reset.
# IDEMPOTENCY GUARD: McGrath reprocess is COMPLETE (mcgrath_persite.parquet analyzed, germline-CpG null). Skip on restart.
sleep 30
[ -f /mnt/data/MCGRATH_WGS_DONE ] && exit 0
GS=gs://ai-temp/apobec-genome-cache; GSUTIL=/usr/bin/gsutil
mkdir -p /mnt/data/scripts
for s in mcgrath_aichem.sh mcgrath_wgs_pileup_aichem.py; do
  for att in 1 2 3 4 5; do $GSUTIL cp $GS/scripts/$s /mnt/data/scripts/$s 2>/dev/null && break; sleep 10; done
done
systemctl is-active mcgrath >/dev/null 2>&1 || systemd-run --unit=mcgrath --property=User=shaharh_quris_ai --setenv=HOME=/home/shaharh_quris_ai /bin/bash -c 'bash /mnt/data/scripts/mcgrath_aichem.sh > /mnt/data/mcgrath.log 2>&1'
