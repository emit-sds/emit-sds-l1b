# David R Thompson
import os

# Intensity-based selection method
if True:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/EMIT_FlatField_20220103'
    cmd = 'python ../utils/makeflat.py --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)

# Manual selection method using a mask image
if False:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/EMIT_FlatField_20220103a'
    cmd = 'python ../utils/makeflat.py --mask_image ../data/emit_radcal_mask.png --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)

