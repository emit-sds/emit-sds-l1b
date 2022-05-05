# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
import os, sys
sys.path.append('../utils')
from fpa import FPA


# This script was used when we were fitting linearity curves to each
# FPA element independently.  It is not used anymore.

#I = envi.open('../data/EMIT_LinearityMap_20220117.hdr').load()

#thresh = 20

#fpa = FPA('../config/tvac2_config.json')


## First version 
#if False:
#    for band in range(I.shape[2]):
#      x = np.squeeze(I[:,:,band])

#      # Remove anomalously high or low values
#      for row in range(1,x.shape[0]):
#        for col in range(x.shape[1]):
#          if abs(x[row,col])>thresh:
#             x[row,col] = x[row-1,col]

#      # Copy and paste linearity columns over the first aquisition zone, 
#      # which is anomalous
#      for col in range(24,44):
#        x[:,col] = x[:,44]

#      # Copy and paste linearity columns over the goober zone, 
#      # which is anomalous
#      for col in range(1020,1027):
#        x[:,col] = x[:,1027]

#      # Copy and paste linearity rows over the OSF filter,
#      # which is anomalous
#      for lo, hi in fpa.osf_seam_positions:
#         for row in range(lo, hi+1):
#           x[row,:] = x[lo-1,:]

#      I[:,:,band] = x.reshape((x.shape[0],x.shape[1],1))

#    envi.save_image('../data/EMIT_LinearityMap_20220117.hdr',I,ext='',force=True)

