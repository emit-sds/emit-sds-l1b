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
import pylab as plt
import json
import logging
import argparse
from numba import jit
from math import pow
from fpa import FPA, frame_embed, frame_extract
from fixosf import get_osf_interp_idx

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')



# Polynomial fitting from https://gist.github.com/kadereub/
@jit(nopython=True)
def _coeff_mat(x, deg):
    mat_ = np.zeros(shape=(x.shape[0],deg + 1))
    const = np.ones_like(x)
    mat_[:,0] = const
    mat_[:, 1] = x
    if deg > 1:
        for n in range(2, deg + 1):
            mat_[:, n] = x**n
    return mat_
    
@jit
def _fit_x(a, b):
    # linalg solves ax = b
    det_ = np.linalg.lstsq(a, b)[0]
    return det_
 
@jit
def fit_poly(x, y, deg):
    a = _coeff_mat(x, deg)
    p = _fit_x(a, y)
    # Reverse order so p[0] is coefficient of highest order
    return p[::-1]

@jit
def eval_polynomial(P, x):
    '''
    Compute polynomial P(x) where P is a vector of coefficients, highest
    order coefficient at P[0].  Uses Horner's Method.
    '''
    result = 0
    for coeff in P:
        result = x * result + coeff
    return result



def spectral_angle(a,b):
    return np.arccos(np.sum(a*b)/(np.sqrt((b**2).sum()) * np.sqrt((a**2).sum())))

def spectral_angle_nocos(a,b):
    return np.sum(a*b)/(np.sqrt((b**2).sum()) * np.sqrt((a**2).sum()))


# Spectral angle comparison
@jit
def closest(a,B):
    numerator = np.sum((B.T*a).T,axis=0)
    denominator = np.sqrt(np.sum((B**2),axis=0)) * np.sqrt(np.sum(a**2))
    projection = numerator/denominator
    return np.argmax(projection)

def ang_infer_bad(frame,t_row,t_col,span,fpa):
    n_rows, n_columns = frame.shape

    closest_col = fpa.SMARGIN
    best_sa = -9e99
    use_idx = np.zeros(n_rows,dtype=bool)
    use_idx[t_row-fpa.BMARGIN:t_row+span+fpa.BMARGIN] = True
    use_idx[t_row:t_row+span] = False

    for col in range(fpa.SMARGIN,n_columns-fpa.SMARGIN):
        if col == t_col:
            continue #Cannot match to itself
        #Get spectral angle

        sa = spectral_angle_nocos(frame[use_idx,t_col],frame[use_idx,col])
        if not np.isnan(sa) and sa > best_sa:
            best_sa = sa
            closest_col = col

    best_spectrum = frame[:,closest_col]
    #slope, offset = fit_poly(best_spectrum[use_idx],frame[use_idx, t_col], 1) 

    sx = best_spectrum[use_idx].sum()
    sy = frame[use_idx, t_col].sum()
    sxy = np.sum(best_spectrum[use_idx]*frame[use_idx, t_col])
    sx2 = np.sum(best_spectrum[use_idx]**2)
    n = np.sum(use_idx)

    slope = ((n*sxy) - (sx*sy))/((n*sx2) - (sx*sx))
    offset = ((sy*sx2) - (sx*sxy))/((n*sx2)-(sx*sx))

    replace_vals = slope * best_spectrum[t_row:t_row+span] + offset
            
    return replace_vals

def ang_replace_bad(frame,bad,fpa):
    n_rows, n_columns = frame.shape
    fixed = frame.copy()

    for col in range(n_columns):
        row = fpa.CMARGIN
        while row < n_rows - fpa.CMARGIN:
            if bad[row,col]>0:
                fixed[row:row+bad[row,col],col] = ang_infer_bad(fixed,row,col,bad[row,col],fpa)
                row += bad[row,col]
            row += 1 #Seems like a bug, but ??      
    return fixed

def ang_fix_bad(frame, bad, fpa):

    rows, columns = frame.shape
    fixed = frame.copy()
    valid_columns = np.where(np.sum(bad,axis=0)==0)[0]
    nfixed = 0

    # Now fix bad pixels in the map using spectral angle 
    # matching plus linear regression inference approach
    # Chapman et al., Remote Sensing 2019
    for col in range(columns):
        if np.sum(bad[:,col])<0:

            #tofix_channels = np.where(bad[:,col]!=0)[0]
            tofix_channels = np.where(bad[:,col]<0)[0] #Make it clear that this code only fixes pixels flagged with <0

            # Don't match on OSF or extrema
            good_channels = (bad[:,col]==0)
            for positions in fpa.osf_seam_positions:
                osf_idx = get_osf_interp_idx(positions)
                good_channels[osf_idx] = False
            good_channels[:fpa.first_illuminated_row] = False
            good_channels[(fpa.last_illuminated_row+1):] = False
            good_channels = np.where(good_channels)[0]

            # calcluate spectral angle over valid FPA elements
            best_sa = 99999
            B = fixed[good_channels, :]
            B = B[:,valid_columns]
            neighbor = closest(fixed[good_channels,col],B)
            best = valid_columns[neighbor]
            best_spectrum = fixed[:,best]
            slope, offset = fit_poly(best_spectrum[good_channels],
                                       fixed[good_channels, col], 1)
            for badc in tofix_channels:
               fixed[badc,col] = slope * best_spectrum[badc] + offset
            nfixed = nfixed + 1
    return fixed


def main():

    description = "Fix bad pixels"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    infile = envi.open(find_header(args.input))
    print(fpa.bad_element_file)
    badfile = envi.open(find_header(fpa.bad_element_file))
    bad = np.squeeze(badfile.load()[:,:,0])

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
            fixed = fix_bad(frame, bad, fpa)

            np.array(fixed, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
