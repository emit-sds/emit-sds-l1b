# David R Thompson

import sys, os
from numpy.linalg import lstsq
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


basedir='/beegfs/scratch/drt/20230901_CPM_GhostStray/20230903_005500_UTC_InFieldStrayLight/'
spatial_files=[basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1100nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1291nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1311nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_1451nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating2_950nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating3_1601nm_darksub_pedestal',
   basedir+'20230903_005500_UTC_InFieldStrayLight_Field180_Grating3_1751nm_darksub_pedestal',
   basedir+'20230903_010000_UTC_InFieldStrayLight_Field180_Grating1_450nm_darksub_pedestal',
   basedir+'20230903_010400_UTC_InFieldStrayLight_Field180_Grating1_600nm_darksub_pedestal',
   basedir+'20230903_010600_UTC_InFieldStrayLight_Field180_Grating1_660nm_darksub_pedestal',
   basedir+'20230903_010900_UTC_InFieldStrayLight_Field180_Grating1_680nm_darksub_pedestal',
   basedir+'20230903_011100_UTC_InFieldStrayLight_Field180_Grating1_800nm_darksub_pedestal',
   basedir+'20230903_013400_UTC_InFieldStrayLight_Field180_Grating3_1901nm_darksub_pedestal']

for infile in spatial_files:
   I = envi.open(infile+'.hdr')
   X = I.load()
   print(X.shape)
   profile = np.median(np.median(X,axis=0),axis=0)
   #plt.plot(profile)
   #plt.show()
   outfile = infile+'_thermcor'
   for i in range(X.shape[0]):
      print(i)
      for j in range(X.shape[1]):
         b = np.squeeze(X[i,j,:])
         coef,_,_,_ = lstsq(profile[:,np.newaxis],b[:,np.newaxis])
         X[i,j,:] = X[i,j,:]-profile*coef
   metadata = I.metadata.copy()
   print(metadata)
   metadata['interleave'] = 'bil'
   envi.save_image(outfile+'.hdr',X,metadata=metadata,interleave='bil',ext='',force=True)
