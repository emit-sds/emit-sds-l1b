# David R Thompson
import os

# Manual selection method using a mask image
if False:
    infile = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_80pc_darksub_pedestal_ghostfix'
    mask = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_80pc_mask.png'
    #backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/CWIS_FlatField_20220218'
    cmd = 'python ../utils/makeflat.py --mask_image %s %s %s' %(mask, infile, outfile)
    print(cmd)
    os.system(cmd)

    infile = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_40pc_darksub_pedestal_ghostfix'
    mask = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_40pc_mask.png'
    #backgroundfile = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_203242_UTC_GenericFOV/20211116_203750_UTC_GenericFOV_Fields-250-1455_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix'
    outfile = '../data/CWIS_FlatField_20220218a'
    cmd = 'python ../utils/makeflat.py --mask_image %s %s %s' %(mask, infile, outfile)
    print(cmd)
    os.system(cmd)

# First radiometric cal
if False:
    infile = '/beegfs/scratch/drt/20220112_CWIS2/20220412_radiometry/20220412_Radiometric_Apo_Box_Scan_darksub_pedestal_ghostfix' 
    outfile = '../data/CWIS_FlatField_20220413'
    cmd = 'python ../utils/makeflat.py --config ../config/cwis2.json %s %s' %(infile, outfile)
    print(cmd)
    os.system(cmd)
 
# Second radiometric cal
for scan in range(1,6):
    infile = '/beegfs/scratch/drt/20220112_CWIS2/20220427_RadCal/20220427_ApoBoxScan_'+str(scan)+'_darksub_pedestal_badfix_osffix_scatterfix_ghostfix'
    outfile = infile+'_flat'#'../data/CWIS_FlatField_20220428'
    cmd = 'python ../../utils/makeflat.py --halfwid 30 --config ../../config/cwis2.json %s %s' %(infile, outfile)
    print(cmd)
    os.system(cmd)
