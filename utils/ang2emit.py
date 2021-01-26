# David R Thompson

import scipy.linalg
import os, sys
import scipy as sp
from spectral.io import envi
import json
import logging
import argparse
import numpy as np
from emit_config import EMITL1Config
from isofit.core.common import resample_spectrum

nchans = 325
nheader = 1

def main():

    description = "Translate coadded AVIRIS-NG radiance data into EMIT format"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('emit_l1_config')
    parser.add_argument('input_rdn')
    parser.add_argument('input_obs')
    parser.add_argument('input_loc')
    parser.add_argument('output_rdn')
    parser.add_argument('output_obs')
    parser.add_argument('output_loc')
    args = parser.parse_args()

    config = EMITL1Config(args.emit_l1_config)

    q, wl, fwhm = np.loadtxt(args.emit_wavelengths).T * 1000.0
    ang_framesize = (640,480)
    I = envi.open(args.input_rdn+'.hdr')
    metadata = I.metadata
    metadata['samples'] = 1280
    metadata['bands'] = nchans + nheader
    wavelength_ang = np.array([float(w) for w in metadata['wavelength']])
    fwhm_ang = np.array([float(f) for f in metadata['fwhm']])
    metadata['wavelength'] = wl 
    metadata['fwhm'] = fwhm
    envi.write_envi_header(output_rdn+'.hdr', metadata)
    
    with open(args.input_rdn,'wb') as finA:
        
        with open(args.input_rdn,'rb') as finA:
        A = sp.fromfile(finA, dtype=sp.float32)

    with open(args.input_B,'rb') as finB:
        B = sp.fromfile(finB, dtype=sp.float32)

    bad = sp.logical_or(A<1e-9,B<1e-9)
    A[bad]=1e-9
    B[bad]=1e-9
    print('Mean absolute difference: %f' % ((abs(A-B)).mean()))

if __name__ == '__main__':

    main()
