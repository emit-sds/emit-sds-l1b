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
from numba import jit


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

nbright = 32

@jit
def addcounts(brightest, frame):
  for row in range(frame.shape[0]):
      for col in range(frame.shape[1]):
          for pos in range(nbright):
             if frame[row,col] > brightest[row,col,pos]:
                 mynext = frame[row,col]
                 # bubble sort
                 for j in range(pos,nbright):
                    swap = brightest[row,col,j]
                    brightest[row,col,j] = mynext
                    mynext = swap
                 break
  return brightest

def main():

    description = "Calculate Flat field"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--cue_channel',default=50,type=int)
    parser.add_argument('--ref_lo',default=99,type=int)
    parser.add_argument('--ref_hi',default=1180,type=int)
    parser.add_argument('--hw_lo',default=8,type=int)
    parser.add_argument('--hw_hi',default=40,type=int)
    parser.add_argument('--selection',type=str,default='spatial')
    parser.add_argument('--badmap_out',type=str,default=None)
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))
 
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
    margin=2

    brightest  = np.ones((rows,columns,nbright))*-9999
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            brightest = addcounts(brightest,frame)
           
        stdev = brightest.std(axis=2)
        flat = brightest.mean(axis=2)
        stdev[np.logical_not(np.isfinite(stdev))] = 0
        flat[np.logical_not(np.isfinite(flat))] = -9999
        reference = np.arange(args.ref_lo,args.ref_hi)
        for row in range(rows):
            ref = np.nanmean(flat[row,reference])
            print('row',row,'reference average is',ref)
            flat[row,:] = flat[row,:] / ref
        flat[np.logical_not(np.isfinite(flat))] = -9999
        envi.save_image(args.output+'.hdr',np.array(flat,dtype=np.float32),ext='',force=True)

if __name__ == '__main__':

    main()
