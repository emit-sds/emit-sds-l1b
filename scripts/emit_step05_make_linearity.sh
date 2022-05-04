# David R Thompson

# collect linearity curves.  Only use the wavelengths from 500-600 nm which are near the photodiode
# this corresponds to a range of channels 274-287
python ../utils/makelinearity.py --plot --top 274 --bottom 287 /beegfs/scratch/drt/20220303_EMIT_TVAC4b/linearity/linearity_Field280*pedestal /beegfs/scratch/drt/20220303_EMIT_TVAC4b/linearity/Field280_clipped_test

# These are vestigial.

#for fieldpoint in 90 165 240 315 390 465 540 615 690 765 840 915 990 1065 1140 1215; do 
#python ../utils/makelinearity.py --top 274 --bottom 287 /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*subframe*pedestal /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_clipped
#done

# collect linearity curves.  Only use the wavelengths from 500-600 nm which are near the photodiode
#for fieldpoint in 290; do 
#python ../utils/makelinearity.py --top 274 --bottom 287 /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*subframe*pedestal /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_clipped
#done

# Collect data points for an example plot at field point 315
#for fieldpoint in 315; do
#python ../utils/makelinearity  --top 274 --bottom 287 /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*subframe*pedestal /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_clipped
#done



