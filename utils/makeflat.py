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
    parser.add_argument('input')
    parser.add_argument('--cue_channel',default=50,type=int)
    parser.add_argument('--ref_lo',default=99,type=int)
    parser.add_argument('--ref_hi',default=1180,type=int)
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

    flat  = np.zeros((rows,columns))
    count = np.zeros((rows,columns))
    allctrs,alllines = [],[]
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            reference = frame[args.cue_channel, :]
            thresh = threshold_otsu(reference)
            lit = np.where(reference>thresh)[0]
            lit.sort()
            ctr = np.median(lit)
            
            # take middle 50% of data
            halfwidth = lit[int(len(lit)/2)]-lit[int(len(lit)/4)]
            if halfwidth<8 or halfwidth>45:
                continue
            for row in range(rows): 
                flat[row,lit] = flat[row,lit] + frame[row,lit] 
                count[row,lit] = count[row,lit] + 1

        flat = flat / count
        reference = np.arange(args.ref_lo,args.ref_hi)
        for row in range(rows):
            flat[row,:] = flat[row,:] / np.mean(flat[row,reference])
        flat[np.logical_not(np.isfinite(flat))] = -9999
        envi.save_image(args.output+'.hdr',np.array(flat,dtype=np.float32),ext='',force=True)

if __name__ == '__main__':

    main()
