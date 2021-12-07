# David R Thompson
import os

#cmd = 'python ../utils/combineflat.py --ref_lo %i --ref_hi %i %s %s' %(infile, ref_lo, ref_hi, outfile)

infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_linear_badfix_scatterfix_ghostfix'
outfile = '../data/EMIT_FlatField_20211207'
ref_lo = 99
ref_hi = 1180
cmd = 'python ../utils/makeflat.py --ref_lo %i --ref_hi %i %s %s' %(ref_lo, ref_hi, infile, outfile)
print(cmd)
os.system(cmd)
