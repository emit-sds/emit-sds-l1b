#! /usr/bin/env python
#
#  Copyright 2021 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import pylab as plt
import json
import logging
import argparse
from numba import jit
from math import pow
from fpa import FPA
from scipy.interpolate import interp1d
import subprocess


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def fix_osf(frame, fpa):
    fixed = frame.copy()

    for positions in fpa.osf_seam_positions:
      osf_idx = get_osf_interp_idx(positions)
      if len(positions) > 2:
        interp_idx = np.array(positions)
        interp_func = interp1d(interp_idx, frame[interp_idx,:], kind='cubic', axis=0)
        osf = np.array(list(set(range(positions[0],positions[-1])).difference(positions)))
      else:
        interp_idx = np.array(positions)
        interp_func = interp1d(interp_idx, frame[interp_idx,:], kind='linear', axis=0)
        osf = np.arange(positions[0]+1, positions[-1])
    fixed[osf_idx,:] = interp_func(osf_idx)

    return fixed

# Determine which indices we're interpolating through for a given OSF seam
def get_osf_interp_idx(positions):
  # If > 2, cubic interpolation
  if len(positions) > 2:
    return np.array(list(set(range(positions[0],positions[-1])).difference(positions)))
  # If = 2, linear interpolation
  else:
    return np.arange(positions[0]+1, positions[-1])



def main():

    description = "Fix bad pixels"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

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

            # Read a frame of data
            if line%10==0:
                logging.info('Line '+str(line))
            print(line)

            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            fixed = fix_osf(frame, fpa)

            np.array(fixed, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
