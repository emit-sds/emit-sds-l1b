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

    description = "Strip a dark frame from an AVIRIS-NG raw file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_raw')
    parser.add_argument('output_dark')
    args = parser.parse_args()

    obc_dark = 2
    rows = 480
    columns = 640
    obc_byte = 641
    dark_margin = 100
    dark_lines = 800
    dark_max = dark_margin + dark_lines
    nframe = rows * columns
    nraw = nframe + columns
    state = 'waiting'
    lines,ndark = 0,0
    darkframes = []

    with open(args.input_raw,'rb') as fin:

        while True:

            # Read a frame of data
            if lines%1000==0:
                logging.info('Calibrating line '+str(lines))
            header = sp.fromfile(fin, count=columns*2, dtype=sp.ubyte)
            frame = sp.fromfile(fin, count=nframe, dtype=sp.uint16)
            frame = frame.reshape((rows, columns))
            
            # Finite state machine
            obc = header[obc_byte]
            if obc == obc_dark:
                ndark = ndark+1

                # If we are in the estimation interval, add this frame
                if ndark > dark_margin and ndark < dark_max:
                    darkframes.append(frame)

                # When we exit the estimation interval, write the average
                # and standard error
                elif ndark >= dark_max:

                    dark_avg = sp.array(darkframes).mean(axis=0)
                    dark_std = sp.array(darkframes).std(axis=0)/sp.sqrt(ndark)

                    with open(args.output_dark,'w') as fout:
                        sp.asarray(dark_avg, dtype=sp.float32).tofile(fout)
                        sp.asarray(dark_std, dtype=sp.float32).tofile(fout)
                    with open(args.output_dark+'.hdr','w') as fout:
                        fout.write(header_string)
                    break

    print('done') 

if __name__ == '__main__':

    main()
