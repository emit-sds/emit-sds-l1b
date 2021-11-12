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
import numpy as np


def main():

    description = "Compare two radiance files"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_A')
    parser.add_argument('input_B')
    args = parser.parse_args()

    with open(args.input_A,'rb') as finA:
        A = sp.fromfile(finA, dtype=sp.float32)

    with open(args.input_B,'rb') as finB:
        B = sp.fromfile(finB, dtype=sp.float32)

    bad = sp.logical_or(A<1e-9,B<1e-9)
    A[bad]=1e-9
    B[bad]=1e-9
    print('Mean absolute difference: %f' % ((abs(A-B)).mean()))

if __name__ == '__main__':

    main()
