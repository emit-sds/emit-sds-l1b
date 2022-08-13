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

    description = 'Multiply two or more flat fields to create another flat';
    parser = argparse.ArgumentParser()

    # Required 
    parser.add_argument('input', nargs='+', help='Flat field')
    parser.add_argument('config')
    parser.add_argument('output', help='Output flat field image')
    args = parser.parse_args()

    fpa = FPA(args.config)

    # Define local variables
    flats = []
    for infile in args.input:
        inhdr  = find_header(infile)
        flat = np.squeeze(envi.open(inhdr).load()[:,:,0])

        if flat.shape[0] < fpa.native_rows or \
           flat.shape[1] < fpa.native_columns:

            # We assume it is clipped
            a = fpa.first_distributed_row
            b = (fpa.last_distributed_row+1)
            c = fpa.first_distributed_column 
            d = (fpa.last_distributed_column+1)
            buffered = np.ones((fpa.native_rows, fpa.native_columns))
            buffered[a:b,c:d] = flat 
            flats.append(buffered)

        else:
            flats.append(flat)

    product = np.prod(np.array(flats),axis=0)
    
    I = envi.open(inhdr)
    meta = I.metadata.copy()
    outhdr = args.output + '.hdr'
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")
    with open(args.output,'wb') as fout:
        np.array(product, dtype=np.float32).tofile(fout)
                   
if __name__ == '__main__':
    main()


