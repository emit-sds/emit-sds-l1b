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
from emitrdn import Config


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

    description = "Synthesize a no-op linearity correction file"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config_file')
    parser.add_argument('output_linearity')
    parser.add_argument('output_bad')
    parser.add_argument('output_rcc')
    parser.add_argument('output_ff')
    parser.add_argument('output_srf')
    parser.add_argument('output_crf')
    parser.add_argument('output_dark')
    args = parser.parse_args()

    with open(config_file,'r') as fin:
        config = json.load(fin)

    for q in [config['radiometric_coefficient_file'],
              config['spectral_calibration_file'],
              config['flat_field_file'],
              config['crf_correction_file'],
              config['srf_correction_file'],
              config['bad_element_file'],
              config['linearity_file'],
              config['dark_frame_file'],
              config['output_rdn_file']]:
        if os.path.exists(q):
             print(q+' exists')
             sys.exit(1)

    # Make Radiometric Coefficients
    n = 238
    rows, cols = 328
    dispersion = 7.4296517
    wl = 265+np.arange(328)*dispersion
    fwhm = (wl[1]-wl[0])*1.1
    rccs = np.array(n) * 0.0001
    uncerts = np.array(n)
    chn = np.arange(n)
    np.savetxt(np.c_[rccs, uncerts, chn], 
        config['radiometric_coefficient_file'])
    np.savetxt(np.c_[chn, wl,fwhm], 
        config['spectral_calibration_file'])
    bad = np.

    with open(args.output_file+'.hdr','w') as fout:



    
    count = sp.arange(2**16, dtype=sp.uint16).tofile(args['linearity_file'])
    with open(config['linearity_file'],'w') as fout:
        fout.write(header_string) 

    print('done') 

if __name__ == '__main__':

    main()
