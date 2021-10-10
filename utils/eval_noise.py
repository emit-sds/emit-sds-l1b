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


header_template = """ENVI
description = {}
samples = 1280
lines = 480
bands = 2
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

rows = 480 
columns = 1280

def main():

    description = "Evaluate noise characteristics"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input', default='')
    parser.add_argument('output', default='')
    parser.add_argument('--start_line', default=0)
    parser.add_argument('--end_line', default=-1)
    args = parser.parse_args()

    in_img = envi.open(args.input+'.hdr')
    
    meta = in_img.metadata.copy()
    nr, nc, nb = meta['lines'],meta['samples'],meta['bands']
    meta['lines'] = rows 
    meta['samples'] = columns 
    meta['bands'] = 2
    meta['data type'] = 4
    meta['interleave'] = 'bsq'

    start_line = args.start_line
    if args.end_line < 0:
        end_line = nr
    out_img = envi.create_file(args.output+'.hdr',metadata=meta,force=True,ext='')
    
    if meta['interleave']=='bil':
        frame = frame.reshape((nr,nb,nc))
    elif meta['interleave']=='bip':
        frame = frame.reshape((nr,nc,nb))
    else:
        raise ValueError('unsupported interleave:',frame)
  
    frames = []
    with open(args.input,'rb') as fin:

        for line in range(nl):
            # Read a frame of data
            if line%1000==0:
                logging.info('Evaluating line '+str(lines))
            header = sp.fromfile(fin, count=columns*2, dtype=sp.ubyte)
            frame = sp.fromfile(fin, count=nframe, dtype=sp.uint16)
            frame = np.array(frame, dtype=np.float32)
            if line >= start_line and line < end_line: 
                frames.append(frame)
    frames= np.array(frames)
    avg = frames.mean(axis=0)       
    stdev=  np.stdev(frames, axis=0)

    with open(args.output, 'wb') as fout:
        np.array(avg,dtype=np.float32).tofile(fout)
        np.array(stdev,dtype=np.float32).tofile(fout)
        
    print('done') 

if __name__ == '__main__':

    main()
