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
from isofit.core.common import svd_inv 
import subprocess


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def precalc_fix_osf_gaussian_variables(fpa, mu, C):

    n_rows = fpa.last_distributed_row - fpa.first_distributed_row + 1

    # Calculate indices for inferred portion and remainder 
    mask = np.ones(n_rows, dtype=bool)
    for positions in fpa.osf_seam_interpolation_edges:
        mask[positions[0]+1:positions[1]] = False
    remain = np.where(mask)[0]
    window = np.where(~mask)[0]

    w = window[:, None]
    r = remain[:, None]
    C11 = C[r, r.T]
    C21 = C[w, r.T]
    Cinv = svd_inv(C11)
    M = C21 @ Cinv

    return {
        "window": window,
        "remain": remain,
        "M": M,
        "mu_window": mu[window],
        "mu_remain": mu[remain]
    }

def fix_osf_gaussian(frame, precomp):
    fixed = frame.copy()
    window = precomp["window"]
    remain = precomp["remain"]
    M = precomp["M"]
    mu_window = precomp["mu_window"]
    mu_remain = precomp["mu_remain"]

    # Perform inference for each spectrum independently
    # mean mu and covariance C are for normalized radiances
    z = np.linalg.norm(frame, axis=0)
    normed = frame / z
    pred = mu_window[:, None] + M @ (normed[remain, :] - mu_remain[:, None])
    fixed[window, :] = pred * z
    return fixed

def fix_osf(frame, fpa):
    fixed = frame.copy()

    if len(fpa.osf_seam_positions)==0:
        return fixed

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

def build_osf_mask(nrows, osf_seam_positions):
    mask = np.ones(nrows, dtype=bool)
    for positions in osf_seam_positions:
        positions_arr = np.asarray(positions, dtype=np.int64)
        idx = get_osf_interp_idx(positions_arr)
        mask[idx] = False
    return mask

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
