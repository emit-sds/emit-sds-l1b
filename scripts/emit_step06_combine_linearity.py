# David R Thompson
import sys, os
import numpy as np
import pylab as plt
from spectral.io import envi


# This is the conservative approach in which the linearity
# correction is a no-op.
if True:

    grid = np.arange(2**16)
    mu = np.ones((1,len(grid)))
    D = np.concatenate((mu,np.zeros((2,len(grid)))),axis=0)

    envi.save_image('../data/EMIT_LinearityBasis_20220504.hdr',
        np.asarray(D,dtype=np.float32),ext='',force=True)

    envi.save_image('../data/EMIT_LinearityMap_20220504.hdr',
        np.zeros((328,1280,2),dtype=np.float32),ext='',force=True)



# This is the expedient approach in which we apply a single lnearity
# curve as a global correction to all FPA elements.
if False:

    curves = []
    for fieldpoint in [280]: 
        fp = str(fieldpoint)
        infile = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/linearity/Field'+fp+'_clipped_test.hdr'
        I = np.squeeze(envi.open(infile).load())
        print(I.shape)
        if len(I.shape)<2:
           I = I.reshape((1,I.size))
        curves.append(I)
    curves = np.concatenate(curves,axis=0)
    grid = np.arange(2**16)
    plt.plot(curves.T)
    p = np.polyfit(grid,curves.mean(axis=0),7)
    mu = np.polyval(p,grid).reshape((1,len(grid)))
    plt.plot(grid,mu.T,'k',linewidth=3)
    plt.show()

    D = np.concatenate((mu,np.zeros((2,len(grid)))),axis=0)

    envi.save_image('../data/EMIT_LinearityBasis_20220423.hdr',
        np.asarray(D,dtype=np.float32),ext='',force=True)

    envi.save_image('../data/EMIT_LinearityMap_20220423.hdr',
        np.zeros((328,1280,2),dtype=np.float32),ext='',force=True)

# This is the old approach in which we create basis vectors
if False:
    # collect linearity curves
    cmd = 'python ../utils/combinelinearity.py '

    # We only use a subset of panel areas to calculate linearity
    #for fieldpoint in [15,90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]:
    for fieldpoint in [165,240,315,390,840,915,990,1065,1140]:
      fp = str(fieldpoint)
      cmd = cmd + '/beegfs/scratch/drt/20211115_EMIT_Linearity/Field_'+fp+'_clipped '
    cmd = cmd + '../data/EMIT_LinearityBasis_20220117'
    print(cmd)
    os.system(cmd)



