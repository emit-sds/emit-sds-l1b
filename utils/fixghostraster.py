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
from fpa import FPA, frame_embed, frame_extract
from fixghost import fix_ghost_matrix
import ray
import pylab as plt


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


@ray.remote
def fix_ghost_parallel(frame, fpa, ghostmap, blur_spatial, blur_spectral):
  return fix_ghost_matrix(frame, fpa, ghostmap, blur_spatial=blur_spatial, 
          blur_spectral=blur_spectral)



# DDA algorithm for line rasterization is
# Courtesy Shivam Pradhan (GeeksForGeeks.org)
def build_ghost_matrix(ghost_config, fpa):

    ghostmap = np.zeros((fpa.native_rows, fpa.native_rows))

    for order in ghost_config['orders']:

      x0, x1 = order['extent']
      scaling = order['scaling']

      # calculate target channel endpoints
      slope, offset = order['slope'], order['offset']
      y0 = x0 * slope + offset 
      y1 = x1 * slope + offset

      # calculate intensity endpoints
      islope, ioffset = order['intensity_slope'], order['intensity_offset']
      i0 = x0 * islope + ioffset
      i1 = x1 * islope + ioffset

      # calculate dx, dy
      dx = x1 - x0;
      dy = y1 - y0;
 
      # calculate steps required for generating pixels
      if abs(dx) > abs(dy): 
          steps = abs(int(round(dx)))
      else:
          steps = abs(int(round(dy)))
 
      # calculate increment in x & y for each step
      xinc = dx / float(steps)
      yinc = dy / float(steps)
 
      # Put pixel for each step
      x = x0
      y = y0
      i = i0
      for j in range(steps+1):
          ghostmap[int(round(x)),int(round(y))] = i
          x = x + xinc
          y = y + yinc
          i = x * islope + ioffset
          i = i * scaling
       
   #plt.imshow(ghostmap)
   #plt.show()
    return ghostmap 


def main():

    description = "Fix spatial and spectral scatter"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--config')
    parser.add_argument('--ncpus',default=30)
    parser.add_argument('ghost_config')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    ray.init()

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

    envi.write_envi_header(args.output+'.hdr',infile.metadata)

    with open(args.ghost_config,'r') as fin:
        ghost_config = json.load(fin)
    ghostmap = build_ghost_matrix(ghost_config, fpa)
    blur_spatial = ghost_config['blur_spatial']
    blur_spectral = ghost_config['blur_spectral']

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

            # embed subframe if needed
            if rows < fpa.native_rows:
                frame = frame_embed(frame, fpa)
            frames.append(frame)

            if len(frames) == args.ncpus or line == (lines-1):
                jobs = [fix_ghost_parallel.remote(f, fpa, ghostmap, blur_spatial,
                                     blur_spectral) for f in frames]
                fixed_all = ray.get(jobs)
                for fixed in fixed_all:

                   # remove embedding if needed
                   if rows < fpa.native_rows:
                       fixed = frame_extract(fixed,fpa)
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
