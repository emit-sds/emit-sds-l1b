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


def noisy_offset(frame, fpa):
    pedestal = np.std(frame[fpa.masked_rows,:],axis=0).mean()
    return pedestal


def main():

    description = "Fix pedestal shift for a data cube"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input', nargs='+')
    parser.add_argument('--config', default=None)
    parser.add_argument('--thresh', default=15)
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    for inpath in args.input:

        infile = envi.open(find_header(inpath))
       
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
        stds = []
       
        with open(inpath,'rb') as fin:
          with open(args.output,'wb') as fout:
       
            for line in range(lines):
       
                # Read a frame of data
                if line%10==0:
                    logging.info('Line '+str(line))
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                stds.append(noisy_offset(frame, fpa))
            stds = np.array(stds) 
            okay = 'OK'
            if np.std(stds)>args.thresh:
               okay = 'BAD'
            outstr = '%s %5.4f %s\n'%(inpath, np.std(stds), okay)
            print(outstr)
            fout.write(bytes(outstr,encoding='utf-8'))
       

if __name__ == '__main__':

    main()
