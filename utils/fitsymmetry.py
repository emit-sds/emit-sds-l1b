# David R Thompson
# Identify the axis of symmetry

import numpy as np
import sys, os
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr
from scipy.optimize import minimize_scalar
from spectral.io import envi

left, right, short, long = 25, 1265, 21, 314
filename = '/beegfs/scratch/drt/20211010_EMIT_Ghost/20211008_114506_FlatField/20211008_114814_UTC_FlatField_Fields-40-1319_darksub'
skip_rows = 400

def symmetric_reflect(img, center):
    rows, cols = img.shape
    new = np.zeros(img.shape)
    for row in range(short,long+1):
        for col in range(left,right+1):
            tcol = center+(center-col)
            neighbors = np.array([np.floor(tcol), np.floor(tcol+1)],dtype=int)
            weights = np.array([1-(tcol-neighbors[0]), 1-(neighbors[1]-tcol)])
            if neighbors[0]>0 and neighbors[0]<1280:
               new[row, neighbors[0]] = new[row, neighbors[0]] + weights[0] * img[row,col]
            if neighbors[1]>0 and neighbors[1]<1280:
               new[row, neighbors[1]] = new[row, neighbors[1]] + weights[1] * img[row,col]
            if row==100 and col==100:
                print(neighbors[0],weights[0])
    return new

def err(x,v):
    r= -np.sum(symmetric_reflect(x,v)*x)
    print(v,':',r)
    return r

ghost = np.zeros((480,480))
X = envi.open(filename+'.hdr').load()

lines, samps, bands = X.shape
half = int(samps/2)
for i in range(skip_rows,lines):
    x = np.squeeze(X[i,:,:]).T
    best = minimize_scalar(lambda v: err(x,v), bounds=[648,655], method='bounded')
    print(best)


