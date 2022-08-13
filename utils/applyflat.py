#!/usr/bin/env python
# David R Thompson 
# Spring 2015
# Jet Propulsion Laboratory, California Institute of Technology

import os, sys, argparse, time
import spectral
from scipy.optimize import minimize
import numpy as np
import spectral.io.envi as envi
from numba import jit
from fpa import FPA, frame_embed, frame_extract


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def main():

    description = 'Apply a multiplicative flat field (and optionally, an additive offset)';
    parser = argparse.ArgumentParser()

    # Required 
    parser.add_argument('input', help='Input radiance image')
    parser.add_argument('config', help='Configuration')
    parser.add_argument('flatfield', help='New (secondary) dark field')
    parser.add_argument('output', help='Output radiance image')
    parser.add_argument('--offset','-r',help='Offset', default=None) 
    args = parser.parse_args()

    # Define local variables
    inhdr  = find_header(args.input)
    I = envi.open(inhdr)
    nrows, ncols, nbands = I.nrows, I.ncols, I.nbands

    outhdr = args.output + '.hdr'
    ffhdr  = args.flatfield+'.hdr'
    flat = np.squeeze(envi.open(ffhdr).load()[:,:,0])

    if args.offset is not None:
        dkhdr  = args.offset+'.hdr'
        offset = np.squeeze(envi.open(dkhdr).load()[:,:,0])
    else:
        offset = np.zeros((nbands, ncols))

    # Clip if necessary
    a,b = fpa.first_distributed_row, (fpa.last_distributed_row+1)
    c,d = fpa.first_distributed_column, (fpa.last_distributed_column+1)
    if flat.shape[0] > nrows:
        flat = flat[a:b,c:d]
    if offset.shape[0] > nrows:
        offset = offset[a:b,c:d]

    meta = I.metadata.copy()
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)
        n = 0

        while frame.size>0:
            frame = frame.reshape((nbands, ncols)) # BIL
            if any((frame<-9989).flatten()):
               continue

            frame = frame * flat + offset 
            frame = np.array(frame,dtype=np.float32)
            frame.tofile(fout)
            frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)
            n = n + 1
            if n%100==0:
                print(n,'/',nrows,' corrected')
                   
if __name__ == '__main__':
    main()


