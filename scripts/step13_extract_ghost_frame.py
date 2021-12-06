# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi

I = envi.open('/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_linear_badfix_scatterfix.hdr')
#I = envi.open('/beegfs/scratch/drt/20211130_EMIT_RadCal/20211116_200400_UTC_GenericFOV/20211116_200836_UTC_GenericFOV_Fields-250-1455_darksub_pedestal_badfix.hdr')
x = I.load()
print(x.shape)
frame = np.squeeze(x[400,:,:])
frame = frame.T
#display = np.log(frame-frame.min()+1e-6)
#display[np.logical_not(np.isfinite(display))] = 1e-6
#plt.semilogy(display[200,:])
#plt.show()
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame.hdr',frame,ext='',force=True)

frame = np.squeeze(x[300,:,:])
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame2.hdr',frame,ext='',force=True)

frame = np.squeeze(x[200,:,:])
frame = frame.T
frame = np.array(frame,dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame3.hdr',frame,ext='',force=True)

print('done')

