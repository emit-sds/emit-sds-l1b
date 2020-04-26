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
description = {Linearity Correction File}
samples = 65536
lines = 1
bands = 0
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""


def main():

    description = "Make a no-op linearity correction file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('output_file')
    args = parser.parse_args()

    count = sp.arange(2**16, dtype=sp.uint16).tofile(args.output_file)
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_string) 

    print('done') 

if __name__ == '__main__':

    main()
