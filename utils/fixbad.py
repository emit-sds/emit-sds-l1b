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


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

def spectral_angle(a,b):
    return np.arccos(np.sum(a*b)/(np.sqrt(pow(a,2).sum()) * np.sqrt(pow(a,2).sum())))


# Spectral angle comparison
def closest(a,B):
    numerator = np.sum((B.T*a).T,axis=0)
    denominator = np.sqrt(np.sum(pow(B,2),axis=0)) * np.sqrt(pow(a,2).sum())
    projection = numerator/denominator
    return np.argmax(projection)


def main():

    description = "Fix bad pixels"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('badmap')
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))
    badfile = envi.open(find_header(args.badmap))
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
            fixed = frame.copy()
            valid_columns = np.where(np.sum(bad,axis=0)==0)[0]
            for col in range(columns):
                if any(bad[:,col]):
                    bad_channels = np.where(bad[:,col]!=0)[0]
                    good_channels = np.where(bad[:,col]==0)[0]
                    best_sa = 99999

                    B = frame[good_channels, :]
                    B = B[:,valid_columns]
                    best = valid_columns[closest(frame[good_channels,col],B)]
                    best_spectrum = frame[:,best]
                   #for test_col in valid_columns:
                   #   sa = spectral_angle(frame[good_channels, test_col],
                   #                       frame[good_channels, col])
                   #   if sa<best_sa:
                   #      best_spectrum = frame[:,test_col]
                   #      best_sa = sa
                    slope, offset = np.polyfit(best_spectrum[good_channels],
                                               frame[good_channels, col], 1)
                    for badc in bad_channels:
                       fixed[badc,col] = slope * best_spectrum[badc] + offset

            np.array(fixed, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
