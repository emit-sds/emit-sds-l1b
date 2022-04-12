# David R Thompson

import sys, os
import numpy as np
from glob import glob
import pylab as plt
from spectral.io import envi


# First compile the smoother library, if needed.
# Go to the ../utils directory and type:
#
# > python setup.py build
#

# Set this to true if we have already completed a full calibration solution
# We can then apply the EMIT calibration process for testing our SRF/CRF 
# correction matrices
validate = True

if False:

    # First combine the data from all the point spread function measurements
    # This requires running basic electronic corrections on all datasets first
    # Record the resulting Gaussian fits in some text files

    files = glob('/beegfs/scratch/drt/20220112_CWIS2/20211210_lasers/20211210_LaserSphere_clip_darksub_pedestal_avg')
    files.sort()
    for target_col in ['84','153','219','276','290','306']:
        for trow in np.arange(40,1280,40):
            target_row = str(trow)
            cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--top_margin','0','--target_col',target_col,'--target_row',target_row]+files+['>>','cwis_spectral_params_clipped.txt']
            print(' '.join(cmds))
            os.system(' '.join(cmds))

if True:

  spatial, spectral = [],[]

  # Now build scatter correction matrices.  Do it twice: first with a null (no-op)
  # correction for comparison, and second for real.
  for magnitude in [0,1]: 

    cmd = 'python ../utils/combinescatter.py --manual '+str(magnitude)+' --spatial  ../scripts/cwis_spectral_params_clipped.txt ../data/CWIS_SpatialScatter_20220406' 
    print(cmd)
    os.system(cmd)

    cmd = 'python ../utils/combinescatter.py --manual '+str(magnitude)+' ../scripts/cwis_spectral_params_clipped.txt ../data/CWIS_SpectralScatter_20220406' 
    print(cmd)
    os.system(cmd)
 
    # Evaluate the result by calibrating a test image
    if validate:

        # Test image for spatial scatter validation 
        dnfile = '/beegfs/scratch/drt/20220112_CWIS2/hiss/20220124t1436avg10_scan500to700nm_clip'
        darkfile = dnfile.replace('_clip','_dark_clip')
        rdnfile = dnfile.replace('_clip','_rdn')
        configfile = '/home/drt/src/emit-sds-l1b/config/cwis2.json'
        cmd = 'python ../emitrdn.py --config %s --dark_file %s %s %s' % (configfile,darkfile,dnfile,rdnfile)
        print(cmd)
        os.system(cmd)
        
        # Extract the point spread function in the spatial dimension
        I = envi.open(rdnfile+'.hdr').load()

        # Plot the result to the screen
        x = np.squeeze(I[69,341,15:])
        spectral.append(x)
        plt.figure(0)
        plt.semilogy(x)
        x = np.squeeze(I[69,:,30])
        spatial.append(x)
        plt.figure(1)
        plt.semilogy(x)

  plt.show()

  # Save pretty plots
  if validate:
      np.savetxt('CWIS_l1bplots_Spatial.txt', np.array(spatial).T)
      np.savetxt('CWIS_l1bplots_Spectral.txt', np.array(spectral).T)
 
 
 
 
 
