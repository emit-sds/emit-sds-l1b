# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from scipy.optimize import minimize
from scipy.interpolate import BSpline,interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def main():

    description = "Calculate Linearity Correction"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('basis')
    parser.add_argument('output')
    args = parser.parse_args()

    xs,ys = [],[]
    nfiles = len(args.input) 
    illums =[] 
    data = []

    basis = envi.open(args.basis+'.hdr').load()
    evec = np.squeeze(basis[1:,:].T)
    if len(evec.shape)<2:
        evec = evec[:,np.newaxis]
    evec[np.isnan(evec)] = 0
    nev = np.squeeze(evec.shape[1])
    mu = np.squeeze(basis[0,:])
    mu[np.isnan(mu)] = 0
    print('mu',mu)
    print('evec',evec)
    data = np.zeros((len(args.input),480,1280))

    for fi,infilepath in enumerate(args.input):

        toks = infilepath.split('_')
        for tok in toks:
            if 'Field' in tok:
               simple = tok.replace('Field','')
               fieldpoint= int(simple)
               active_cols = np.arange(fieldpoint-37,fieldpoint+38,dtype=int)
            elif 'candelam2' in tok:
               simple = tok.split('.')[0]
               simple = simple.replace('PD','')
               simple = simple.replace('candelam2','')
               simple = simple.replace('p','.')
               illums.append(float(simple))
        
        infile = envi.open(find_header(infilepath))
        
        if int(infile.metadata['data type']) == 2:
            dtype = np.uint16
        elif int(infile.metadata['data type']) == 4:
            dtype = np.float32
        else:
            raise ValueError('Unsupported data type')
        if infile.metadata['interleave'] != 'bil':
            raise ValueError('Unsupported interleave')
        
        rows = int(infile.metadata['bands'])
        columns = int(infile.metadata['samples'])
        lines = int(infile.metadata['lines'])
        nframe = rows * columns
        
        sequence = []
        with open(infilepath,'rb') as fin:
        
            for line in range(lines):
        
                # Read a frame of data
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                sequence.append(frame[:,active_cols])
                
        sequence = np.array(sequence)
        data[fi,:,active_cols] = np.median(sequence, axis=0).T
               
    out = np.zeros((480,1280,nev))
    for wl in np.arange(26,313):
   
       for col in range(columns):
       
           DN = data[:,wl,col]
           L = np.array(illums) 
           L = L / L.mean() * DN.mean()
           
           # best least-squares slope forcing zero intercept
           tofit = np.where(np.logical_and(DN>1000,DN<35000))[0]
           slope = np.sum(DN[tofit]*L[tofit])/np.sum(pow(L[tofit],2))
           if not np.isfinite(slope):
              continue
           ideal = slope * L
   
           # Don't correct above the saturation level
           ideal[DN>40000] = DN[DN>40000]
           
           grid = np.arange(2**16)
           resamp = interp1d(DN, ideal, bounds_error=False, fill_value='extrapolate')(grid)
           resamp = resamp / grid
           resamp[np.isnan(resamp)]=0
           coef = (resamp - mu) @ evec 
         # if wl>30 and col>100 and col<1200:
         #     plt.plot(resamp)
         #     plt.plot(np.squeeze(evec@coef[:,np.newaxis]) + mu,'k.')
         #     plt.show()
         # out[wl,col,:] = coef
         # print('!',wl,col,coef)

    envi.save_image(args.output+'.hdr',np.array(out,dtype=np.float32),ext='',force=True)

if __name__ == '__main__':

    main()
