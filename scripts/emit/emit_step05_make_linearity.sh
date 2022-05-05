# David R Thompson

# collect linearity curves.  Only use the wavelengths from 500-600 nm which are near the photodiode
# this corresponds to a range of channels 274-287
python ../utils/makelinearity.py --plot --top 274 --bottom 287 /beegfs/scratch/drt/20220303_EMIT_TVAC4b/linearity/linearity_Field280*pedestal /beegfs/scratch/drt/20220303_EMIT_TVAC4b/linearity/Field280_clipped_test


