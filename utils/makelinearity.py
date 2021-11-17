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
    parser.add_argument('--cue_channel',type=int,default=150)
    args = parser.parse_args()

    xs,ys = [],[]
    
    for infilepath in args.input:

        illum = float(infilepath.split('_')[-3][2:-9].replace('p','.'))
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
                DN = reference[1215]
                x.append(DN)
                y.append(illum)
               
            ys.append(np.median(y))
            xs.append(np.median(x))

    plt.plot(ys,xs,'ko')
    
    plt.xlabel('DN')
    plt.ylabel('candela m2')
    plt.box(False)
    plt.grid(True)
    plt.ylim([0,5000])
    plt.xlim([0,500])

    use = xs<40000
    p = np.polyfit(xs[use],ys[use],1)
    plt.figure()
    plt.plot(xs[use],np.polyval(p,xs[use])-ys[use],'ko')
    plt.box(False)
    plt.grid(True)
    plt.show()

if __name__ == '__main__':

    main()
