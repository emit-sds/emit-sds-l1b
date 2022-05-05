# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all 
# the SRF data from the monochromator sweeps
# It outputs text files, one per base file, corresponding to the SRFs

# Make bad pixels
for fin in `ls /beegfs/scratch/drt/20211113_EMIT_SRF/all/*subframe_darksub_pedestal`; do
   python ../../utils/makesrf.py $fin > ${fin}.txt
done


