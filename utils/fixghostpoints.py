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
import ray


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

@ray.remote
def fix_ghost(frame, ghostmap, center=649.5, blur_spatial=50, blur_spectral=1, fudge = 4.25):

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
                 ghost[ghost_row, tcol] + frame[row,col] * intensity * fudge

  start = 25
  ghost[start:,:] = gaussian_filter(ghost[start:,:],[blur_spectral, blur_spatial])
  new = frame - ghost
  return new



def main():

    description = "Fix spatial and spectral scatter"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--ncpus',default=30)
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

    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    lines = int(infile.metadata['lines'])
    nframe = rows * columns

    ray.init()

    envi.write_envi_header(args.output+'.hdr',infile.metadata)

    ghost_config = np.loadtxt(args.ghost_config)
    ghostmap = [[] for r in range(rows)]
    for source, target, confidence, intensity in ghost_config:
        ghostmap[int(source)].append((int(target),intensity))

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        frames = []
        for line in range(lines):

            # Read a frame of data
            if line%10==0:
                print('Line '+str(line))

            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            if infile.metadata['interleave'] == 'bil':
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            elif infile.metadata['interleave'] == 'bip':
                frame = np.array(frame.reshape((columns,rows)),dtype=np.float32).T
            else:
                raise ValueError('unsupported interleave')

            frames.append(frame)

            if len(frames) == args.ncpus or line == (lines-1):
                jobs = [fix_ghost.remote(f, ghostmap) for f in frames]
                fixed_all = ray.get(jobs)
                for fixed in fixed_all:
                   if infile.metadata['interleave'] == 'bil':
                       np.array(fixed, dtype=np.float32).tofile(fout)
                   elif infile.metadata['interleave'] == 'bip':
                       np.array(fixed.T, dtype=np.float32).tofile(fout)
                   else:
                       raise ValueError('unsupported interleave')
                frames = []

    print('done') 

if __name__ == '__main__':

    main()