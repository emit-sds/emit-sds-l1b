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
    parser.add_argument('input', nargs='+', help='Flat field')
    parser.add_argument('output', help='Output radiance image')
    args = parser.parse_args()

    # Define local variables
    flats = []
    for infile in args.input:
        inhdr  = find_header(infile)
        flat = np.squeeze(envi.open(inhdr).load()[:,:,0])
        flats.append(flat)
    avgflat = np.median(np.array(flats),axis=0)
    
    I = envi.open(inhdr)
    meta = I.metadata.copy()
    outhdr = args.output + '.hdr'
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")
    with open(args.output,'wb') as fout:
        np.array(avgflat, dtype=np.float32).tofile(fout)
                   
if __name__ == '__main__':
    main()


