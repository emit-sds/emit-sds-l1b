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


left, right, top, bottom = 25, 1264, 26, 313

def main():

    description = "Calculate Flat field"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('output')
    args = parser.parse_args()

    xs,ys = [],[]
    nfiles = len(args.input) 
    illums =[] 
    data = []

    for fi,infilepath in enumerate(args.input):

        toks = infilepath.split('_')
        for tok in toks:
            if 'Field' in tok:
               simple = tok.replace('Field','')
               fieldpoint= int(simple)
               active_rows = np.arange(max(left,fieldpoint-37),min(fieldpoint+38,right),dtype=int)
            elif 'candelam2' in tok:
               simple = tok.split('.')[0]
               simple = simple.replace('PD','')
               simple = simple.replace('candelam2','')
               simple = simple.replace('p','.')
               illums.append(float(simple))
        print(infilepath)
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
        
        x,y = [],[]
        with open(infilepath,'rb') as fin:
        
            image_data = np.zeros((rows,columns))
            for line in range(lines):
        
                # Read a frame of data
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                image_data[:,active_rows] = frame[:,active_rows]
        data.append(image_data)        

    data = np.array(data)
    curves = []
    for wl in np.arange(top,bottom):
        DN = data[:,wl,active_rows].mean(axis=1)
        L = np.array(illums) 

        L = L / L.mean() * DN.mean()
        
        # best least-squares slope forcing zero intercept
        tofit = np.where(np.logical_and(DN>1000,DN<35000))[0]
        use = np.where(np.logical_and(DN>10,DN<35000))[0]

        slope = np.sum(DN[tofit]*L[tofit])/np.sum(pow(L[tofit],2))
        print(slope)   
        if not np.isfinite(slope):
           continue
        ideal = slope*L

        # Don't correct above the saturation level
        ideal[DN>40000] = DN[DN>40000]
        
        smoothed = lowess(ideal,DN,frac=0.4,return_sorted=False)
        smoothed[DN>1000]=DN[DN>1000]

        grid = np.arange(2**16)
        resamp = interp1d(DN, ideal, bounds_error=False, fill_value='extrapolate')(grid)
        curves.append(resamp)

    curves = np.array(curves,dtype=np.float32)
    envi.save_image(args.output+'.hdr',curves,ext='',force=True)

    if False:
            plt.figure(0)
            #plt.loglog(grid,curves[np.arange(0,curves.shape[0],50),:].T,'k-')
            plt.plot(grid,(curves[np.arange(0,curves.shape[0],50),:]/grid).T,'k-')
            #plt.semilogx(ys[use],ideal[use]/corrected[use],color+'-')
            #plt.plot(xscaled[use],ideal[use],color+'+')
            #plt.semilogx(xs[use],(corrected[use]-target[use])/(target[use])*100,'ro')
            #plt.semilogx(ys[use],np.zeros(len(use)),'k:')
           #plt.xlabel('DN')
           #plt.ylabel('Deviation from linearity (%)')
            plt.box(False)
            plt.grid(True)

           #plt.figure(1)
           #plt.semilogx(grid,resamp/grid,'k')
            #plt.plot(xscaled[use],ideal[use],color+'+')
            #plt.semilogx(xs[use],(corrected[use]-target[use])/(target[use])*100,'ro')
            #plt.semilogx(ys[use],np.zeros(len(use)),'k:')
           #plt.xlabel('DN')
           #plt.ylabel('Deviation from linearity (%)')
           #plt.box(False)
           #plt.grid(True)

            plt.show()

if __name__ == '__main__':

    main()
