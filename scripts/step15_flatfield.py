# David R Thompson
import os

if False:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix'
    outfile = '../data/EMIT_FlatField_20211216'
    cmd = 'python ../utils/makeflat.py --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)

if True:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_linear'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_linear'
    outfile = '../data/EMIT_FlatField_20211228'
    cmd = 'python ../utils/makeflat.py --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)
