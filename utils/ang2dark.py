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
description = {AVIRIS-NG Dark Frame}
samples = 640
lines = 480
bands = 2
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""


def main():

    description = "Strip a dark frame from an AVIRIS-NG file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_file', nargs=1, default='')
    parser.add_argument('output_file', nargs=1, default='')
    args = parser.parse_args()

    OBC_DARK = 2
    rows = 480
    columns = 640
    nframe = rows * columns
    nraw = nframe + colummns
    state = 'waiting'
    lines,ndark = 0,0
    with open(args.input_file,'rb') as fin:

        # Read a frame of data
        if lines%1000==0:
            logging.info('Calibrating line '+str(lines))
        header = sp.fromfile(fin, count=columns*2, dtype=sp.ubyte)
        frame = sp.fromfile(fin, count=nframe, dtype=sp.uint16)
        
        # Finite state machine
        obc = header[obc_byte]
        if ndark < dark_margin:
            if obc == OBC_DARK:
                ndark = ndark+1
            if ndark > 100 and ndark <= 800:
                darkframes.append(frame)
            elif ndark>900:
                dark_avg = s.array(darkframes).mean(axis=0)
                dark_std = s.array(darkframes).std(axis=0)/s.sqrt(ndark)
                with open(args.input_file,'rb') as fin:
                    sp.asarray(dark_avg, dtype=s.float32).tofile(fout)
                    sp.asarray(dark_std, dtype=s.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
