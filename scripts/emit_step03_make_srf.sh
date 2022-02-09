# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all the SRF data

# Make bad pixels
for fin in `ls /beegfs/scratch/drt/20211113_EMIT_SRF/all/*clip*pedestal`; do
   python ../utils/makesrf.py $fin > ${fin}.txt
done


