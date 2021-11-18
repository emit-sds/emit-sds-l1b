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
from scipy.interpolate import splrep,splev
from skimage.filters import threshold_otsu
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def main():

    description = "Calculate Flat field"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('--cue_channel',type=int,default=270)
    args = parser.parse_args()

    xs,ys = [],[]
    
    for infilepath in args.input:

        toks = infilepath.split('_')
        for tok in toks:
            if 'Field' in tok:
               simple = tok.replace('Field','')
               fieldpoint= int(simple)
            elif 'candelam2' in tok:
               simple = tok.split('.')[0]
               simple = simple.replace('PD','')
               simple = simple.replace('candelam2','')
               simple = simple.replace('p','.')
               illum = float(simple)
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
        
            for line in range(lines):
        
                # Read a frame of data
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                reference = frame[args.cue_channel, :]
                DN = reference[fieldpoint+1]
                x.append(DN)
                y.append(illum)
               
            xs.append(np.median(x))
            ys.append(np.median(y))

    plt.plot(ys,xs,'ko')
    
    plt.xlabel('DN')
    plt.ylabel('candela m2')
    plt.box(False)
    plt.grid(True)

    xs = np.array(xs,dtype=float)
    ys = np.array(ys,dtype=float)
    target = ys / ys.mean() * xs.mean()
    use = np.where(np.logical_and(xs>5,xs<40000))[0]
    tofit = np.where(np.logical_and(xs>5,xs<40000))[0]

    # best least-squares slope forcing zero intercept
    slope = np.sum(ys[tofit]*xs[tofit])/np.sum(pow(xs[tofit],2))
    slope = slope/ys[tofit].mean() * xs[tofit].mean()
    print(slope)
  
    # Third order correction, see Fiedler et al., 
    # APPLIED OPTICS  Vol. 44, No. 25  1 September 2005

    def model(x, DN):
      maxdn = 40000
      return x[0]*DN/maxdn + x[1]*pow(DN/maxdn,2) + x[2]*pow(DN/maxdn,3)

    def err(x, DN, tru):
      return np.mean(pow((model(x,DN) - tru)/tru,2))

   #x0 = np.array([40000,0,0])
   #best = minimize(lambda x: err(x,xs,target), x0)
    p = splrep(xs, target, k=3, s=4)
    print(p)

    if True:
        #corrected = model(best.x, xs)
        corrected = splev(xs,p)
        plt.figure()
        plt.semilogx(xs[use],(xs[use]-target[use])/(target[use])*100,'ko')
        plt.semilogx(xs[use],(corrected[use]-target[use])/(target[use])*100,'ro')
        plt.semilogx(xs[use],np.zeros(len(use)),'k:')
        plt.xlabel('DN')
        plt.ylabel('Deviation from linearity (%)')
        plt.box(False)
        plt.grid(True)
        plt.show()

if __name__ == '__main__':

    main()
