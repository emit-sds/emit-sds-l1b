# David R Thompson

import sys, os
import numpy as np
from glob import glob
import pylab as plt
from spectral.io import envi

config = '../../config/cpm.json'

# First compile the smoother library, if needed.
# Go to the ../utils directory and type:
#
# > python setup.py build
#

if True:
  basedir = '/beegfs/scratch/drt/20230901_CPM_GhostStray/20230901_060455_UTC_SPECTRALSTRAYLIGHT/'
  spectral_files = [\
                basedir+'20230901_060630_UTC_SpectralStrayLight_Field180_Grating1_450nm_darksub_pedestal', 
                basedir+'20230901_060806_UTC_SpectralStrayLight_Field180_Grating1_600nm_darksub_pedestal',
                #basedir+'20230901_060912_UTC_SpectralStrayLight_Field180_Grating1_660nm_darksub_pedestal',
                #basedir+'20230901_061017_UTC_SpectralStrayLight_Field180_Grating1_680nm_darksub_pedestal',
                basedir+'20230901_061122_UTC_SpectralStrayLight_Field180_Grating1_800nm_darksub_pedestal',
                #basedir+'20230901_061227_UTC_SpectralStrayLight_Field180_Grating2_950nm_darksub_pedestal',
                basedir+'20230901_061352_UTC_SpectralStrayLight_Field180_Grating2_1101nm_darksub_pedestal',
                basedir+'20230901_061456_UTC_SpectralStrayLight_Field180_Grating2_1292nm_darksub_pedestal',
                basedir+'20230901_061602_UTC_SpectralStrayLight_Field180_Grating2_1311nm_darksub_pedestal',
                basedir+'20230901_061707_UTC_SpectralStrayLight_Field180_Grating2_1451nm_darksub_pedestal',
                basedir+'20230901_061813_UTC_SpectralStrayLight_Field180_Grating3_1601nm_darksub_pedestal',
                basedir+'20230901_061935_UTC_SpectralStrayLight_Field180_Grating3_1751nm_darksub_pedestal',
                basedir+'20230901_062043_UTC_SpectralStrayLight_Field180_Grating3_1901nm_darksub_pedestal',
                basedir+'20230901_062148_UTC_SpectralStrayLight_Field180_Grating3_2051nm_darksub_pedestal',
                basedir+'20230901_062253_UTC_SpectralStrayLight_Field180_Grating3_2200nm_darksub_pedestal',
                basedir+'20230901_062357_UTC_SpectralStrayLight_Field180_Grating3_2349nm_darksub_pedestal',
                basedir+'20230901_062500_UTC_SpectralStrayLight_Field180_Grating3_2498nm_darksub_pedestal']

  # ,'--plot',
  cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--target_row','180']+\
       spectral_files+['>','../../data/cpm/spectral_params_clipped.txt']
  print(' '.join(cmds))
  os.system(' '.join(cmds))

  # Now build scatter correction matrices.  Do it twice: first with a null (no-op)
  # correction for comparison, and second for real.
  magnitude = 1
  cmd = 'python ../../utils/combinescatter.py --constant  --config '+config+' --manual '+str(magnitude)+\
        ' ../../data/cpm/spectral_params_clipped.txt ../../data/cpm/CPM_SpectralSccatter_20231002' 
  print('---------------')
  print(cmd)
  os.system(cmd)

 
if False:
  basedir='/beegfs/scratch/drt/20230901_CPM_GhostStray/20230903_005500_UTC_InFieldStrayLight/'
  spatial_files=[\
     basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1100nm_darksub_pedestal_thermcor',
     basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1291nm_darksub_pedestal_thermcor',
     basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1311nm_darksub_pedestal_thermcor',
     #basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1451nm_darksub_pedestal_thermcor',
     basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_950nm_darksub_pedestal_thermcor',
     #basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating3_1601nm_darksub_pedestal_thermcor',
     #basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating3_1751nm_darksub_pedestal_thermcor',
     #basedir+'20230903_010000_UTC_InFieldStrayLight_Field180_Grating1_450nm_darksub_pedestal_thermcor',
     #basedir+'20230903_010400_UTC_InFieldStrayLight_Field180_Grating1_600nm_darksub_pedestal_thermcor',
     #basedir+'20230903_010600_UTC_InFieldStrayLight_Field180_Grating1_660nm_darksub_pedestal_thermcor',
     #basedir+'20230903_010900_UTC_InFieldStrayLight_Field180_Grating1_680nm_darksub_pedestal_thermcor',
     basedir+'20230903_011100_UTC_InFieldStrayLight_Field180_Grating1_800nm_darksub_pedestal_thermcor']
     #basedir+'20230903_013400_UTC_InFieldStrayLight_Field180_Grating3_1901nm_darksub_pedestal_thermcor']

  spatial_files.sort()
  #,'--plot'
  cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--min_signal','4e-5','--target_row','180','--spatial']+\
         spatial_files+['>','../../data/cpm/spatial_params_clipped.txt']
  print(' '.join(cmds))
  os.system(' '.join(cmds))
  
  magnitude = 1
  # Use constant relationship w/ wavelength
  cmd = 'python ../../utils/combinescatter.py --constant --config '+config+' --manual '+str(magnitude)+\
        ' --spatial  ../../data/cpm/spatial_params_clipped.txt ../../data/cpm/CPM_SpatialScatter_20231002' 
  print('---------------')
  print(cmd)
  os.system(cmd)



 
 
 
 
