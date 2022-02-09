#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import json
import logging
import argparse
from numba import jit
from math import pow
from fpa import FPA, frame_embed, frame_extract



def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def fix_scatter(frame, spectral_correction, spatial_correction):
   if frame.shape[0] != spectral_correction.shape[0] or \
       frame.shape[1] != spatial_correction.shape[1]:
       logging.error('Mismatched frame size')
   fixed = spectral_correction @ (spatial_correction @ frame.T).T
   return fixed


def main():

    description = "Fix spatial and spectral scatter"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--config')
    parser.add_argument('spatial_corr')
    parser.add_argument('spectral_corr')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    infile = envi.open(find_header(args.input))
    spatialfile = envi.open(find_header(args.spatial_corr))
    spatial = np.squeeze(spatialfile.load())
    spectralfile = envi.open(find_header(args.spectral_corr))
    spectral = np.squeeze(spectralfile.load())

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

            # Read a frame of data
            if line%10==0:
                logging.info('Line '+str(line))

            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            fixed = fix_scatter(frame, spectral, spatial)
            np.array(fixed, dtype=np.float32).tofile(fout)

 
    print('done') 

if __name__ == '__main__':

    main()
