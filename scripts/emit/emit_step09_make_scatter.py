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

    files = glob('/beegfs/scratch/drt/20211114_EMIT_Infield/20211114_InFieldScatter/*clip*linear')
    files.sort()
    cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--target_row','40','--spatial']+files+['>>','spatial_params_clipped.txt']
    os.system(' '.join(cmds))


    files = glob('/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/*clip*linear')
    files.sort()
    cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--target_row','940']+files+['>>','spectral_params_clipped.txt']
    os.system(' '.join(cmds))

if True:

  spatial, spectral = [],[]

  # Now build scatter correction matrices.  Do it twice: first with a null (no-op)
  # correction for comparison, and second for real.
  for magnitude in [0,1]: 

    cmd = 'python ../../utils/combinescatter.py --manual '+str(magnitude)+' --spatial  ../../scripts/spatial_params_clipped.txt ../../data/EMIT_SpatialScatter_20220406' 
    print(cmd)
    os.system(cmd)

    cmd = 'python ../../utils/combinescatter.py --manual '+str(magnitude)+' ../../scripts/spectral_params_clipped.txt ../../data/EMIT_SpectralScatter_20220406' 
    print(cmd)
    os.system(cmd)
 
    # Evaluate the result by calibrating a test image
    if validate:

        # Test image for spatial scatter validation 
        testdir = '/beegfs/scratch/drt/20211114_EMIT_Infield/20211114_InFieldScatter/'
        darkfile = testdir + '20211114_051100_UTC_InFieldScatter_dark_clip'
        dnfile = testdir + '20211114_051117_UTC_InFieldScatter_2058p14nm_clip'
        config = '../../config/EMIT_20220504.json'
        rdnfile = dnfile + '_rdn' #dnfile.replace('.raw','_rdn')
        cmd = 'python ../../emitrdn.py %s %s %s %s' % (dnfile,config,darkfile,rdnfile)
        print(cmd)
        os.system(cmd)
        
        # Extract the point spread function in the spatial dimension
        I = envi.open(rdnfile+'.hdr').load()
        I = np.squeeze(np.mean(I,axis=0)).T
        band = np.argmax(np.mean(I,axis=1))
        I = np.squeeze(I[band,:])
        
        # Plot the result to the screen
        plt.figure(0)
        spatial.append(I)
        plt.semilogy(I)

        # Test image for spectral scatter validation 
        testdir = '/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/'
        darkfile = testdir+'20211115_225554_UTC_InFieldScatter_dark_clip'
        dnfile = testdir+'20211115_225605_UTC_InFieldScatter_2060p12nm_clip'
        rdnfile = dnfile + '_rdn'#.replace('.raw','_rdn')
        cmd = 'python ../emitrdn.py --dark_file %s %s %s' % (darkfile,dnfile,rdnfile)
        print(cmd)
        os.system(cmd)
 
        # Extract the point spread function in the spectral dimension
        I = envi.open(rdnfile+'.hdr').load()
        I = np.squeeze(np.mean(I,axis=0)).T
        position = np.argmax(np.mean(I,axis=0))
        I = np.squeeze(I[:,position])
        
        # Plot the result to the screen
        plt.figure(1)
        spectral.append(I)
        plt.semilogy(I)

  plt.show()

  # Save pretty plots
  if validate:
      np.savetxt('../data/plots/EMIT_l1bplots_Spatial.txt', np.array(spatial).T)
      np.savetxt('../data/plots/EMIT_l1bplots_Spectral.txt', np.array(spectral).T)
 
 
 
 
 
