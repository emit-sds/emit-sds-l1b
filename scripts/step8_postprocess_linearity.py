# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi

I = envi.open('../data/EMIT_LinearityMap_20211215.hdr').load()

thresh = 20


for band in range(I.shape[2]):
  x = np.squeeze(I[:,:,band])

  # Remove anomalously high or low values
  for row in range(1,x.shape[0]):
    for col in range(x.shape[1]):
      if abs(x[row,col])>thresh:
         x[row,col] = x[row-1,col]

  # Copy and paste linearity columns over the first aquisition zone, 
  # which is anomalous
  for col in range(24,44):
    x[:,col] = x[:,44]

  # Copy and paste linearity columns over the goober zone, 
  # which is anomalous
  for col in range(1020,1027):
    x[:,col] = x[:,1027]

  # Copy and paste linearity rows over the OSF filter,
  # which is anomalous
  for row in range(187,189):
    x[row,:] = x[186,:]

  I[:,:,band] = x.reshape((x.shape[0],x.shape[1],1))

envi.save_image('../data/EMIT_LinearityMap_20211215.hdr',I,ext='',force=True)


# collect new linearity curves
#for fieldpoint in 15 90 165 240 315 390 465 540 615 690 765 840 915 990 1065 1140 1215; do 
#for fieldpoint in 690 765 840 915 990 1065 1140 1215; do 

# python ~/src/emit-sds-l1b/utils/makelinearity.py --plot /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*linear /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_linearcheck

#done


# Matador test!

#python ../utils/matador.py \
#/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_025005_UTC_LinearitySphere_Field1215_Step14p8mm_PD1289p0candelam2_darksub_pedestal_badfix \
#/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_024003_UTC_LinearitySphere_Field1215_Step18p4mm_PD6p0candelam2_darksub_pedestal_badfix \
#/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_025005_UTC_LinearitySphere_Field1215_Step14p8mm_PD1289p0candelam2_darksub_pedestal_badfix_linear \
#/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_024003_UTC_LinearitySphere_Field1215_Step18p4mm_PD6p0candelam2_darksub_pedestal_badfix_linear
