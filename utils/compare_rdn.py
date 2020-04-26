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



def main():

    description = "Compare two files"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_A')
    parser.add_argument('input_B')
    args = parser.parse_args()

    with open(args.input_A,'rb') as finA:
        A = sp.fromfile(finA, dtype=sp.float32)
        print(A)

    with open(args.input_B,'rb') as finB:
        B = sp.fromfile(finB, dtype=sp.float32)
        print(B)

    print('Mean absolute difference: %f' % (abs(A-B).mean()))

if __name__ == '__main__':

    main()
