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
from scipy.ndimage import gaussian_filter
from numba import jit
from math import pow


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

@jit   
def fix_ghost(frame, ghostmap, center=649.5, blur_spatial=3, blur_spectral=3):

  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  if rows>cols:
      raise IndexError('Misformed frame')

  for row in range(rows):
    for ghost_row, intensity in ghostmap[row]:
       for col in range(cols):
          tcol = int(center*2 - col)
          if tcol>0 and tcol<1280:
              ghost[ghost_row, tcol] = \
                 ghost[ghost_row, tcol] + frame[row,col] * intensity

  ghost = gaussian_filter(ghost,[blur_spectral, blur_spatial])
  new = frame - ghost
  return new



def main():

    description = "Fix spatial and spectral scatter"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('ghost_config')
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))

    if int(infile.metadata['data type']) == 2:
        dtype = np.uint16
    elif int(infile.metadata['data type']) == 4:
        dtype = np.float32
    else:
        raise ValueError('Unsupported data type')

    with open(args.ghost_config,'r') as fin:
        config = json.load(fin)

    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    lines = int(infile.metadata['lines'])
    nframe = rows * columns

    envi.write_envi_header(args.output+'.hdr',infile.metadata)

    ghost_config = np.loadtxt(args.ghost_config)
    ghostmap = [[] for r in rows]
    for source, target, confidence, intensity in ghost_config:
        ghostmap[source].append(target,intensity)

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        for line in range(lines):

            # Read a frame of data
            if line%10==0:
                print('Line '+str(line))

            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            if infile.metadata['interleave'] == 'bil':
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                fixed = fix_ghost(frame, ghostmap)
                np.array(fixed, dtype=np.float32).tofile(fout)
            elif infile.metadata['interleave'] == 'bip':
                frame = np.array(frame.reshape((columns,rows)),dtype=np.float32).T
                fixed = fix_ghost(frame, ghostmap)
                np.array(fixed.T, dtype=np.float32).tofile(fout)
            else:
                raise ValueError('unsupported interleave')

    print('done') 

if __name__ == '__main__':

    main()
