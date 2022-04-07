# David R Thompson
import os, sys, glob
from spectral.io import envi
import numpy as np

# This creates a no-op linearity basis and coefficient map
grid = np.arange(2**16)
basis = np.concatenate((np.ones((1,len(grid))),
                        np.zeros((1,len(grid))),
                        np.zeros((1,len(grid)))),axis=0)
envi.save_image('../data/CWIS_LinearityBasis_22020331.hdr',
    np.array(basis[:,:,np.newaxis],dtype=np.float32),
    ext='',force=True,interleave='bsq')
 
cmap = np.zeros((328,1280,2))
envi.save_image('../data/CWIS_LinearityMap_22020331.hdr',
    np.array(cmap, dtype=np.float32),
    ext='',force=True,interleave='bip')

