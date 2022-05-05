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
from fpa import FPA, frame_embed, frame_extract

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def fix_pedestal(frame, fpa):
    if fpa.pedestal_strategy == 'neither':
       total, counts = [],[]
       mask = np.zeros(frame.shape)>1
       if len(fpa.masked_cols)>0:
          mask[:,fpa.masked_cols] = True
       if len(fpa.masked_rows)>0:
          mask[fpa.masked_rows,:] = True
       avg = frame[mask].mean() * fpa.pedestal_multiplier 
       frame = frame - avg
       return frame
    if (fpa.pedestal_strategy == 'columns' or \
          fpa.pedestal_strategy == 'both') and \
          len(fpa.masked_cols)>0:
        pedestal = np.median(frame[:,fpa.masked_cols], axis=1)
        pedestal = pedestal * fpa.pedestal_multiplier
        frame = (frame.T-pedestal).T
    if (fpa.pedestal_strategy == 'rows' or \
          fpa.pedestal_strategy == 'both') and \
          len(fpa.masked_rows)>0:
        pedestal = np.median(frame[fpa.masked_rows,:], axis=0)
        pedestal = pedestal * fpa.pedestal_multiplier
        frame = frame-pedestal
    return frame


def main():

    description = "Fix pedestal shift for a data cube"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    infile = envi.open(find_header(args.input))

    if int(infile.metadata['data type']) == 2:
        dtype = np.int16
    elif int(infile.metadata['data type']) == 12:
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


    metadata = infile.metadata.copy()
    metadata['data type'] = 4
    envi.write_envi_header(args.output+'.hdr', metadata)

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        for line in range(lines):

            # Read a frame of data
            if line%10==0:
                logging.info('Line '+str(line))
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            fixed = fix_pedestal(frame, fpa)
            np.array(fixed, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
