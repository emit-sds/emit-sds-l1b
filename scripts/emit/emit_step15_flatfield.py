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
if False:
   infile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal'
   outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_flat'
   backgroundfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal'
   cmd = 'python ../utils/makeflat.py --halfwid 30 --background %s %s %s' %(backgroundfile, infile, outfile)
   print(cmd)
   os.system(cmd)

# TVAC4B sphere sequence at OAP focus
if False:
   infile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/asd_full_infocus/emit20220308t200138_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal'
   outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/asd_full_infocus/emit20220308t200138_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal_flat'
   #outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/flatfield_infocus/emit20220307t205822_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal_flat'
   cmd = 'python ../utils/makeflat.py --halfwid 30 %s %s' %(infile, outfile)
   print(cmd)
   os.system(cmd)

# TVAC4B sphere sequence at OAP focus
if False:
   infile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/asd_half_infocus/emit20220309t003507_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal'
   outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/asd_half_infocus/emit20220309t003507_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal_flat'
   #outfile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/flatfield_infocus/emit20220307t205822_o00000_s000_l1a_raw_b0101_v01_strip_shift_darksub_pedestal_flat'
   cmd = 'python ../utils/makeflat.py --halfwid 30 %s %s' %(infile, outfile)
   print(cmd)
   os.system(cmd)

infiles = ['/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_50/emit20220305t022009_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_50/emit20220305t022009_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_50_half/emit20220305t191035_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_50_half/emit20220305t191035_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_half/emit20220305t161821_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_half/emit20220305t161821_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix']

unused =  ['/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_dark/emit20220305t033519_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_dark/emit20220305t033519_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix']

bgfiles = ['/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked/emit20220305t012240_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked_half/emit20220305t174320_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked_half/emit20220305t174320_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked_half/emit20220305t174320_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix',
           '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/radcal_blocked_half/emit20220305t174320_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix']

for infile, bgfile in zip(infiles, bgfiles):

   outfile = infile + '_flat'
   cmd = 'srun -n 1 -N 1 -c 40 --mem 180GB python ../utils/makeflat.py --background %s --halfwid 30 %s %s' %(bgfile, infile, outfile)
   print(cmd)
   #os.system(cmd)


basedir= '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal_nonuniform/'

cmd = 'cp ' 
cmd = cmd+basedir+'/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix_flat'
cmd = cmd + ' ' + '../data/EMIT_FlatFieldHalf_20220504'
print(cmd)
os.system(cmd)

cmd = 'cp ' 
cmd = cmd+basedir+'/radcal/emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix_flat.hdr'
cmd = cmd + ' ' + '../data/EMIT_FlatFieldHalf_20220504.hdr'
print(cmd)
os.system(cmd)

cmd = 'cp ' 
cmd = cmd+basedir+'/radcal_half/emit20220305t161821_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix_flat'
cmd = cmd + ' ' + '../data/EMIT_FlatField_20220504'
print(cmd)
os.system(cmd)

cmd = 'cp ' 
cmd = cmd+basedir+'/radcal_half/emit20220305t161821_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_badfix_osffix_scatterfix_ghostfix_flat.hdr'
cmd = cmd + ' ' + '../data/EMIT_FlatField_20220504.hdr'
print(cmd)
os.system(cmd)


