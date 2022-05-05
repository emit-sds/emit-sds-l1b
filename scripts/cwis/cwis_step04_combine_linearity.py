# David R Thompson
import sys, os
import numpy as np
from spectral.io import envi

# collect linearity curves
cmd = 'python ~/src/emit-sds-l1b/utils/combinelinearity.py /beegfs/scratch/drt/20220112_CWIS2/linearity/cwis_linearity_curves ../data/CWIS2_LinearityBasis_20220115'
print(cmd)
#os.system(cmd)

# Now fit per wavelength
cmd = 'python ~/src/emit-sds-l1b/utils/fitlinearity.py --config ../../config/cwis2.json --width 10 --margin 0 /beegfs/scratch/drt/20220112_CWIS2/linearity/PD_Field365_PD*candelam2_darksub_pedestal ../../data/CWIS2_LinearityBasis_20220115 ../../data/CWIS2_LinearityMap_20220116'
print(cmd)
#os.system(cmd)

if True:
    I = envi.open('../../data/CWIS2_LinearityMap_20220116.hdr').load()

    for band in range(I.shape[2]):
      x = np.squeeze(I[:, 365, band])
      print(x.shape)
      for i in range(I.shape[1]):
         I[:,i,band] = x.reshape((len(x),1,1))#np.tile(x, (1, I.shape[1], 1))

    envi.save_image('../../data/CWIS2_LinearityMap_20220116.hdr', I, ext='', force=True)


