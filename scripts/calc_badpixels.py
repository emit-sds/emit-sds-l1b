# David R Thompson

from spectral.io import envi
import numpy as np
import os,sys
import pylab as plt

'''
This script analyzes a standard dark data cube to identify bad pixels
of the EMIT FPA.
'''

dark_file = '20210817_161811_ARF/20210817_161820_UTC_ARF_dark'
base = os.getenv('EMIT_TVAC_DATA_DIRECTORY')
X = envi.open(base+'/'+dark_file+'.hdr')
X = X.load()
use_lines = 50
X = X[:use_lines,:,:]

mu = np.squeeze(X.mean(axis=0)).T
nc,ns = mu.shape

mu = mu.reshape((nc*ns))
thresh = 1000
bad = abs(mu-np.median(mu))>thresh
bad = bad.reshape(nc,ns)

bads = 0
bad_map = bad.copy()
bad_map = np.array(bad_map,dtype=np.int16)
for column in range(bad_map.shape[1]):
    state_machine = 0
    for row in range(bad_map.shape[0]):
        if bad[row,column]:
            state_machine = state_machine + 1
            bad_map[row,column] = -state_machine
            print(row,column,state_machine)
            if row<328:
                bads = bads + 1
        else:
            state_machine = 0
print('total bads:',bads)
bad_map = bad_map.reshape((nc,ns,1))
envi.save_image('../data/EMIT_Bad_Elements_20210817.hdr',
        bad_map, interleave='bsq', ext='', force=True)

