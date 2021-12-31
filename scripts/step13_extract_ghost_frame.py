# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi


raddir = '/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/'
I = envi.open(raddir+'20211115_195942_UTC_InFieldScatter_1003p64nm_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr')
x = I.load()
frame = np.squeeze(x.mean(axis=0))
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame.hdr',frame,ext='',force=True)

raddir = '/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/'
I = envi.open(raddir+'20211115_204013_UTC_InFieldScatter_1249p16nm_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr')
x = I.load()
frame = np.squeeze(x.mean(axis=0))
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame2.hdr',frame,ext='',force=True)

lindir = '/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/'
I = envi.open(lindir + '20211117_084042_UTC_LinearitySphere_Field240_Step15p2mm_PD1000p0candelam2_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr')
x = I.load()
frame = np.squeeze(x.mean(axis=0))
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame3.hdr',frame,ext='',force=True)

lindir = '/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/'
I = envi.open(lindir + '20211117_044932_UTC_LinearitySphere_Field840_Step14p4mm_PD1616p0candelam2_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr')
x = I.load()
frame = np.squeeze(x.mean(axis=0))
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame4.hdr',frame,ext='',force=True)

raddir = '/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/'
I = envi.open(raddir+'20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr')
x = I.load()
frame = np.squeeze(x[340:350,:,:].mean(axis=0))
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame5.hdr',frame,ext='',force=True)
print('done')

