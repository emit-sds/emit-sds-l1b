# David R Thompson
# david.r.thompson@jpl.nasa.gov

# Define filepaths.  We should already have dark-subtracted and pedestal-shift-corrected the input
basedir=/beegfs/scratch/drt/20220303_EMIT_TVAC4b/flatfield/
original=${basedir}/emit20220304t173343_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal

# First average the flat field scene
averaged=${basedir}/emit20220304t173343_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_avg
python ../../utils/emit2dark.py ${original} ${averaged}

# Apply the bad pixel mask to the average 
output=../data/EMIT_BadElements_20220307
python ../../utils/makebad.py ${averaged} ${output}
