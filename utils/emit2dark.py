#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import scipy as sp
from spectral.io import envi
import json
import logging
import argparse


header_string = """ENVI
description = {EMIT Dark Frame}
lines = %i
samples = %i
bands = 2
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def dark_from_file(filepath):

    infile = envi.open(find_header(filepath))

    if int(infile.metadata['data type']) == 2:
        dtype = sp.int16
    elif int(infile.metadata['data type']) == 12:
        dtype = sp.uint16
    elif int(infile.metadata['data type']) == 4:
        dtype = sp.float32
    else:
        raise ValueError('Unsupported data type')
    if infile.metadata['interleave'] != 'bil':
        raise ValueError('Unsupported interleave')

    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    lines = int(infile.metadata['lines'])
    nframe = rows * columns

    darkframes = []
    with open(filepath, 'rb') as fin:
        for line in range(lines):
        
            # Read a frame of data
            if line%10==0:
                logging.info('Averaging line '+str(line))
            frame = sp.fromfile(fin, count=nframe, dtype=dtype)
            if len(frame)<nframe:
                print('exiting early')
                break
            frame = frame.reshape((rows, columns))
            darkframes.append(frame)

    ndark = len(darkframes)
    dark_avg = sp.array(darkframes).mean(axis=0)

    if ndark>1:
        dark_std = sp.array(darkframes).std(axis=0)/sp.sqrt(ndark)
    else:
        dark_std = np.ones((dark_avg.shape)) * -9999
    return dark_avg, dark_std


  

def main():

    description = "Average a dark sequence with outlier rejection"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('output_dark')
    args = parser.parse_args()

    dark_avg, dark_std = dark_from_file(args.input)

    with open(args.output_dark,'w') as fout:
        sp.asarray(dark_avg, dtype=sp.float32).tofile(fout)
        sp.asarray(dark_std, dtype=sp.float32).tofile(fout)

    with open(args.output_dark+'.hdr','w') as fout:
        fout.write(header_string % dark_avg.shape)

    print('done') 

if __name__ == '__main__':

    main()
