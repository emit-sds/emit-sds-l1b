# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from numba import jit


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


@jit
def fix_linearity(frame, mu, evec, coeffs):
    rows, columns = frame.shape
    new = np.zeros(frame.shape) 
    for row in range(rows):
        for col in range(columns):
            i = max(int(frame[row,col]),0)
            tot = np.sum(evec[i,:]*coeffs[row,col,:]) + mu[i]
            new[row,col] = frame[row,col] * tot
    return new


def main():

    description = "Linearity correction"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('basis')
    parser.add_argument('coefficients')
    parser.add_argument('output')
    args = parser.parse_args()

    xs,ys = [],[]
    nfiles = len(args.input) 
    illums =[] 
    data = []

    basis = envi.open(args.basis+'.hdr').load()
    evec = np.squeeze(basis[1:,:].T)
    evec[np.isnan(evec)] = 0
    mu = np.squeeze(basis[0,:])
    mu[np.isnan(mu)] = 0

    coeffs = envi.open(args.coefficients+'.hdr').load()
   
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

    envi.write_envi_header(args.output+'.hdr',infile.metadata)
    
    with open(args.input,'rb') as fin:
        with open(args.output,'wb') as fout:
        
            for line in range(lines):
        
                print(line)
                # Read a frame of data
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                new = fix_linearity(frame, mu, evec, coeffs)
                np.array(new,dtype=np.float32).tofile(fout)

if __name__ == '__main__':

    main()
