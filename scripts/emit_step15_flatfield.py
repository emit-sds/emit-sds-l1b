# David R Thompson
import os

# Intensity-based selection method for TVAC2 data
if False:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/EMIT_FlatField_20220117'
    cmd = 'python ../utils/makeflat.py --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)

# Manual selection method using a mask image for TVAC2 data
if False:
    infile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/EMIT_FlatField_20220117'
    cmd = 'python ../utils/makeflat.py --mask_image ../data/emit_radcal_mask.png --background %s %s %s' %(backgroundfile, infile, outfile)
    print(cmd)
    os.system(cmd)

# TVAC4B panel sequence
if True:
   infile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal'
   outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_flat'
   backgroundfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal'
   cmd = 'python ../utils/makeflat.py --background %s %s %s' %(backgroundfile, infile, outfile)
   print(cmd)
   os.system(cmd)


