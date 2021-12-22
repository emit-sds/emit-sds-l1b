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

    description = "Make linearity curves"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--chan_lo',type=int)
    parser.add_argument('--chan_hi',type=int)
    parser.add_argument('--fullrange',action='store_true')
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

    # Set wavength ranges used for the estimation
    lo, hi = top, bottom
    if args.chan_lo is not None:
        lo = args.chan_lo
    if args.chan_hi is not None:
        hi = args.chan_hi
        
    for wl in np.arange(lo,hi):

        DN = data[:,wl,active_rows].mean(axis=1)
        L = np.array(illums) 

        L = L / L.mean() * DN.mean()
        
        # best least-squares slope forcing zero intercept
        tofit = np.where(np.logical_and(DN>1000,DN<35000))[0]
        tofit = np.where(np.logical_and(DN>100,DN<30000))[0]

        slope = np.sum(DN[tofit]*L[tofit])/np.sum(pow(L[tofit],2))
        #slope, offset = np.polyfit(L[tofit],DN[tofit],1)
        print(slope)#,offset) 

        if not np.isfinite(slope):
           continue

        ideal = slope*L
        grid = np.arange(2**16)

        resamp = interp1d(DN, ideal, bounds_error=False, fill_value='extrapolate')(grid)
        resamp = resamp / grid

        if not args.fullrange:
            resamp[grid<1000]=resamp[np.argmin(abs(grid-1000))]
        # Don't correct above the saturation level
        resamp[grid>40000]=resamp[np.argmin(abs(grid-40000))]

       #checkband = np.argmin(abs(grid-20000))
       #if resamp[checkband]>0.998 and resamp[checkband]<1.002:
        if True:
            curves.append(resamp)
    curves = np.array(curves,dtype=np.float32)
    print(curves.shape)
    envi.save_image(args.output+'.hdr',curves,ext='',force=True)

    if args.plot:
            plt.figure(0)
            plt.semilogx(DN, (DN/ideal).T,'k.')
            #plt.semilogx(grid,(curves[np.arange(0,curves.shape[0],50),:]).T,'k-')
            plt.box(False)
            plt.grid(True)
            plt.ylim([0.9,3.0])
            plt.show()

if __name__ == '__main__':

    main()
