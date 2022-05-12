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
from math import pow
from fpa import FPA
from fixghost import fix_ghost
import ray
import pylab as plt
from scipy.stats import norm


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


@ray.remote
def fix_ghost_parallel(frame, fpa, ghostmap, blur, center, plot):
  return fix_ghost(frame, fpa, ghostmap, blur=blur, center=center, plot=plot)



# DDA algorithm for line rasterization is
# Courtesy Shivam Pradhan (GeeksForGeeks.org)
def build_ghost_matrix(ghost_config, fpa):

    ghostmap = np.zeros((fpa.native_rows, fpa.native_rows))

    for order in ghost_config['orders']:

      x0, x1 = order['extent']
      x0 = max(x0,fpa.first_valid_row)
      x1 = min(x1,fpa.last_valid_row)
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
          if int(round(x))<=fpa.last_valid_row and\
             int(round(y))<=fpa.last_valid_row:
              ghostmap[int(round(x)),int(round(y))] = i
          x = x + xinc
          y = y + yinc
          i = x * islope + ioffset
          i = i * scaling
       
   #plt.imshow(ghostmap)
   #plt.show()
    return ghostmap 


def build_ghost_blur(ghost_config, fpa):

    blur_dictionary = {}
    for area in ghost_config['psf_zones']:

        indices = np.arange(area['extent'][0],area['extent'][1]+1)
        spatial_blur = np.zeros((fpa.native_columns, fpa.native_columns))
        spectral_blur = np.zeros((fpa.native_rows, fpa.native_rows))

        for psf in area['psfs']:

            sigma = psf['sigma']
            peak = psf['peak']

            # Spectral blur over specific range
            for chan in indices:
                lineshape = norm.pdf(indices,chan,sigma)
                lineshape = lineshape/max(lineshape) * peak
                spectral_blur[chan,indices] = spectral_blur[chan,indices] + lineshape 

            # Spatial blur over entire range
            cols = range(fpa.native_columns)
            for col in cols:
                lineshape = norm.pdf(cols,col,sigma)
                lineshape = lineshape/max(lineshape) * peak
                spatial_blur[col,:] = spatial_blur[col,:] + lineshape 

        # Add the identity
        spectral_blur = spectral_blur + np.eye(fpa.native_rows)
        spatial_blur = spatial_blur + np.eye(fpa.native_columns)

        # Should we normalize? Scaling is redundant with 
        # Ghost intensity parameter.  
        for chan in range(spectral_blur.shape[0]):
             spectral_blur[chan,:] = spectral_blur[chan,:]/sum(spectral_blur[chan,:])
        for col in range(spatial_blur.shape[0]):
             spatial_blur[col,:] = spatial_blur[col,:]/sum(spatial_blur[col,:])
                
        blur_dictionary[(area['extent'][0], area['extent'][1])] = \
                [spatial_blur, spectral_blur]

    return blur_dictionary


def main():

    description = "Fix spatial and spectral scatter"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--ncpus',default=30)
    parser.add_argument('config')
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

    with open(fpa.ghost_map_file,'r') as fin:
        ghost_config = json.load(fin)
    ghostmap = build_ghost_matrix(ghost_config, fpa)
    blur = build_ghost_blur(ghost_config, fpa)
    center = ghost_config['center']

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
            frames.append(frame)

            if len(frames) == args.ncpus or line == (lines-1):
                jobs = [fix_ghost_parallel.remote(f, fpa, ghostmap, 
                                     blur, center=center, plot=args.plot) for f in frames]
                fixed_all = ray.get(jobs)
                for fixed in fixed_all:

                   # remove embedding if needed
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
