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

    description = 'Combine multiple flat fields with a median';
    parser = argparse.ArgumentParser()

    # Required 
    parser.add_argument('input', type=str, help='Text file with flat field paths')
    parser.add_argument('output', help='Output radiance image')
    args = parser.parse_args()

    # Define local variables
    flats = []
    with open(args.input,'r') as fin:
        for line in fin.readlines():
            inhdr  = find_header(line.strip())
            flat = np.squeeze(envi.open(inhdr).load()[:,:,0])
            flats.append(flat)
    median_flat = np.median(np.array(flats),axis=0)
    
    I = envi.open(inhdr)
    meta = I.metadata.copy()
    outhdr = args.output + '.hdr'
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")
    with open(args.output,'wb') as fout:
        np.array(median_flat, dtype=np.float32).tofile(fout)
                   
if __name__ == '__main__':
    main()


