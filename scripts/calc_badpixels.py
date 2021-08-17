# David R Thompson

from spectral.io import envi
import numpy as np
import os,sys
import pylab as plt

dark_file = '20210817_155151_CRF/20210817_155200_UTC_CRF_dark'
base = os.getenv('EMIT_TVAC_DATA_DIRECTORY')
X = envi.open(base+'/'+dark_file+'.hdr')
X = X.load()
use_lines = 50
X = X[:use_lines,:,:]

mu = np.squeeze(X.mean(axis=0)).T
nc,ns = mu.shape

mu = mu.reshape((nc*ns))
bad = abs(mu-np.median(mu))>1000
bad = bad.reshape(nc,ns)


bad_map = bad.copy()
bad_map = np.array(bad_map,dtype=np.int16)
for column in range(bad_map.shape[1]):
    state_machine = 0
    for row in range(bad_map.shape[0]):
        if bad[row,column]:
            state_machine = state_machine + 1
            bad_map[row,column] = -state_machine
            print(row,column,state_machine)
        else:
            state_machine = 0
bad_map = bad_map.reshape((nc,ns,1))
envi.save_image('../data/EMIT_Bad_Elements_20210817.hdr',
        bad_map, interleave='bsq', ext='', force=True)

