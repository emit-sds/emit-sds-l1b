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
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
from makelinearity import linearize
from fpa import FPA
import scipy.linalg as linalg
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
    parser.add_argument('--config')
    parser.add_argument('--linearity_nbasis',default=2)
    parser.add_argument('--width',default=37)
    parser.add_argument('--margin',default=9)
    parser.add_argument('--draft',default=None)
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    margin = int(args.margin)
    width = int(args.width)
    xs,ys = [],[]
    nfiles = len(args.input) 
    illums =[] 
    out = np.zeros((fpa.native_rows,fpa.native_columns,args.linearity_nbasis))
    if args.draft is not None:
        out = envi.open(args.draft+'.hdr').load()

    basis = np.squeeze(envi.open(args.basis+'.hdr').load())
    evec = np.squeeze(basis[1:,:].T)
    if evec.shape[1] != args.linearity_nbasis:
        raise IndexError('Linearity basis does not match file size')
    evec[np.isnan(evec)] = 0
    for i in range(args.linearity_nbasis):
      evec[:,i] = evec[:,i] / linalg.norm(evec[:,i])
    print(linalg.norm(evec,axis=1),linalg.norm(evec,axis=0))
    mu = np.squeeze(basis[0,:])
    mu[np.isnan(mu)] = 0
    data, last_fieldpoint = [], -9999


    for fi,infilepath in enumerate(args.input):

        print('loading %i/%i: %s'%(fi,len(args.input),infilepath))
        toks = infilepath.split('_')
        for tok in toks:
            if 'Field' in tok:
               simple = tok.replace('Field','')
               fieldpoint= int(simple)
               if last_fieldpoint<0: 
                   last_fieldpoint = fieldpoint
               elif last_fieldpoint != fieldpoint:
                   raise IndexError('One fieldpoint per call. Use --draft')
               margin=9
               active_cols = np.arange(max(fpa.first_illuminated_column, fieldpoint-width-margin),
                                       min(fpa.last_illuminated_column ,fieldpoint+width+1-margin),dtype=int)
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

        infile = envi.open(infilepath+'.hdr')
        frame_data = infile.load().mean(axis=0)
        data.append(frame_data[active_cols,:])
    data = np.array(data) 
               
    for wl in np.arange(fpa.first_valid_row, fpa.last_valid_row+1):
   
       for mycol,col in enumerate(active_cols):

         DN = data[:,mycol,wl]
         L = np.array(illums) 
         resamp = linearize(DN, L )#,plot=(wl>50 and col>40 and col<1200))
         coef = (resamp - mu)[np.newaxis,:] @ evec
         out[wl,col,:] = coef[:args.linearity_nbasis]
         if False:#wl>50 and col>40 and col<1200:
             plt.plot(resamp)
             plt.plot(resamp-mu)
             plt.plot(np.squeeze(np.sum(evec*coef,axis=1)) + mu,'k.')
             plt.show()
         if wl%10==0:
             print('!',wl,col,coef)

    envi.save_image(args.output+'.hdr',np.array(out,dtype=np.float32),ext='',force=True)

if __name__ == '__main__':

    main()
