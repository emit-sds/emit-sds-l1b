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
    parser.add_argument('--cue_channel',default=148,type=int)
    parser.add_argument('--ref_lo',default=99,type=int)
    parser.add_argument('--ref_hi',default=1180,type=int)
    parser.add_argument('--hw_lo',default=50,type=int)
    parser.add_argument('--hw_hi',default=180,type=int)
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

    flat  = np.zeros((rows,columns))
    count = np.zeros((rows,columns))
    sumsq = np.zeros((rows,columns))
    ref = np.zeros((lines,columns))
    allctrs,alllines = [],[]
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            ref[line,:] = frame[args.cue_channel, :]

    thresh = np.sort(ref,axis=0)
    thresh = thresh[-10,:]

    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            reference = frame[args.cue_channel, :]
            use = np.where(reference>thresh)[0]

            print(line,len(use),np.median(use),thresh)
            flat[:,use] = flat[:,use] + frame[:,use] 
            count[:,use] = count[:,use] + 1
            sumsq[:,use] = sumsq[:,use] + pow(frame[:,use],2)

        mean_sumsq = sumsq / count
        flat = flat / count

        rowmean = flat[:,30:1250].mean(axis=1)
        rowstdev = flat[:,30:1250].std(axis=1)
        stdev = np.sqrt(mean_sumsq - pow(flat,2))
        stdev[np.logical_not(np.isfinite(stdev))] = 0
        bad = np.logical_or(np.logical_or(stdev==0,
              (abs(flat.T-rowmean)>rowstdev*20).T),stdev>100)
    
    bad[:,:25] = 0
    bad[:,1265:] = 0
   #plt.hist(stdev.flatten(),500)
   #plt.figure()
   #plt.imshow(bad)
   #plt.show()
    bads = 0
    bad_map = bad.copy()
    bad_map = np.array(bad_map,dtype=np.int16)
    for column in range(bad_map.shape[1]):
        state_machine = 0
        for row in range(bad_map.shape[0]):
            if bad[row,column]:
                state_machine = state_machine + 1
                bad_map[row,column] = -state_machine
                print(row,column,state_machine)
                if row<328:
                    bads = bads + 1
            else:
                state_machine = 0
    print('total bads:',bads)
    bad_map = bad_map.reshape((rows,columns,1))
    envi.save_image(args.output+'.hdr',
        bad_map, interleave='bsq', ext='', force=True)

if __name__ == '__main__':

    main()
