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
from fpa import FPA
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


# Create linearity curve for a paired list of DNs and illuminations
def linearize(DN, L, plot=False):

     L = L / L.mean() * DN.mean()
          
     # best least-squares slope forcing zero intercept
     tofit = np.where(np.logical_and(DN>2000,DN<36000))[0]
     if len(tofit)<1:
        return np.ones((2**16),dtype=float)
     slope = np.sum(DN[tofit]*L[tofit])/np.sum(pow(L[tofit],2))

     ideal = slope*L
     grid = np.arange(2**16)

     # First we extrapolate the linearity behavior at the bottom 5% of the
     # dynamic range, overwriting measured DNs
     extrap_range = DN<2000
     extrap_train = np.logical_and(DN>2000,DN<8000)
     if sum(extrap_train)>1 and sum(extrap_range)>0:
         p = np.polyfit(DN[extrap_train], 
                        DN[extrap_train]/ideal[extrap_train], 1)
         DN[extrap_range] = np.polyval(p, DN[extrap_range]) * ideal[extrap_range]
      
     extrap_range = DN>42000
     extrap_train = np.logical_and(DN>30000,DN<43000)
     if sum(extrap_train)>1 and sum(extrap_range)>0:
         p = np.polyfit(DN[extrap_train], 
                        DN[extrap_train]/ideal[extrap_train], 1)
         DN[extrap_range] = np.polyval(p, DN[extrap_range]) * ideal[extrap_range]
     else:
         print(DN)
     resamp = interp1d(DN, ideal, bounds_error=False, fill_value='extrapolate')(grid)
     resamp = resamp / grid
     resamp[np.logical_not(np.isfinite(resamp))] = 1.0
     resamp[:25] = resamp[25] 

     if plot:
        plt.figure(0)
        plt.plot(DN, 1.0/(DN/ideal).T,'k.')
        plt.plot(grid, resamp,'r')
        plt.box(False)
        plt.grid(True)
        plt.ylim([0.9,3.0])
        np.savetxt('linearity_plot_points.txt', np.c_[DN,1.0/(DN/ideal).T])
        np.savetxt('linearity_plot_interp.txt', np.c_[grid,resamp])
        plt.show()


     return resamp



def main():

    description = "Make linearity curves"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--margin',default=37)
    parser.add_argument('--config')
    parser.add_argument('output')
    args = parser.parse_args()
    margin = int(args.margin)

    fpa = FPA(args.config)
    left = fpa.first_illuminated_column
    right = fpa.last_illuminated_column
    top = fpa.first_illuminated_row
    bottom = fpa.last_illuminated_row

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
               active_rows = np.arange(max(left,fieldpoint-margin),
                                      min(fieldpoint+margin,right)+1,dtype=int)
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
        
        x,y = [],[]
        image_data = (infile.load())[:,active_rows,:].mean(axis=1).mean(axis=0)
        print(image_data.shape)
        data.append(image_data)

    data = np.array(data)
    curves = []
    print(data.shape)
        
    for wl in np.arange(top,bottom):

        DN = data[:,wl]
        L = np.array(illums) 
        resamp = linearize(DN, L, plot=(args.plot and wl==100))
        if all(np.logical_and(resamp>0.98,resamp<1.02)):
            curves.append(resamp)
  
    curves = np.array(curves,dtype=np.float32)
    envi.save_image(args.output+'.hdr',curves,ext='',force=True)

    
if __name__ == '__main__':

    main()
