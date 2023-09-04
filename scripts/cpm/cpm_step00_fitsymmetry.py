# David R Thompson
# Identify the axis of symmetry

import numpy as np
import sys, os
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr
from scipy.optimize import minimize_scalar, brute
from spectral.io import envi

left, right, short, long = 20, 620, 40, 400
filename = '/beegfs/scratch/drt/20230901_CPM_GhostStray/20230901_072600_UTC_GhostMapping_Whitelight/20230901_073900_UTC_GhostMapping_Whitelight_OG590_darksub_pedestal'
skip_rows = 0

def symmetric_reflect(img, center):
    rows, cols = img.shape
    new = np.zeros(img.shape)
    for col in range(left,right+1):
        tcol = int(center+(center-col))
        if tcol>=0 and tcol<640:
            new[:, tcol] = img[:,col]
    return new

def err(x,v):
    r= -np.sum(symmetric_reflect(x,v)*x)
    #print(v,':',r)
    return r

ghost = np.zeros((480,480))
X = envi.open(filename+'.hdr').load()

lines, samps, bands = X.shape
half = int(samps/2)
for i in range(skip_rows,lines):
    x = np.squeeze(X[i,:,:]).T
    best = brute(lambda v: err(x,v), ranges=(slice(300,340,1),))
    print(best)


