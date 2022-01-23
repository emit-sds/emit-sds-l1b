#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys, os.path
import scipy as sp
import numpy as np
from spectral.io import envi
from datetime import datetime, timezone
from scipy import linalg, polyfit, polyval
import json
import logging
import argparse
import multiprocessing
import ray

# Import some EMIT-specific functions
my_directory, my_executable = os.path.split(os.path.abspath(__file__))
sys.path.append(my_directory + '/utils/')

from fpa import FPA, frame_embed, frame_extract
from fixbad import fix_bad
from fixosf import fix_osf
from fixlinearity import fix_linearity
from fixscatter import fix_scatter
from fixghost import fix_ghost
from fixghostraster import build_ghost_matrix
from pedestal import fix_pedestal
from darksubtract import subtract_dark
from emit2dark import dark_from_file


header_template = """ENVI
description = {{Calibrated Radiance, microWatts per (steradian nanometer [centemeter squared])}}
samples = {ncolumns}
lines = {lines}
bands = {nchannels}
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
byte order = 0
wavelength units = Nanometers
wavelength = {{{wavelength_string}}}
fwhm = {{{fwhm_string}}}
band names = {{{band_names_string}}}"""

 
def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

     
class Config:

    def __init__(self, fpa, filename, dark_file=None):

        # Load calibration file data
        with open(filename,'r') as fin:
         self.__dict__ = json.load(fin)
        
        # Adjust local filepaths where needed
        for fi in ['spectral_calibration_file','srf_correction_file',
                   'crf_correction_file','linearity_file','ghost_map_file',
                   'radiometric_coefficient_file', 'linearity_map_file',
                   'bad_element_file','flat_field_file']:
            path = getattr(self,fi)
            if path[0] != '/':
                path = os.path.join(my_directory, path)
                setattr(self,fi,path)

        if dark_file is not None:
            self.dark_frame_file = dark_file
        self.dark, self.dark_std = dark_from_file(self.dark_frame_file)
        _, self.wl_full, self.fwhm_full = \
             sp.loadtxt(self.spectral_calibration_file).T * 1000
        self.srf_correction = sp.fromfile(self.srf_correction_file,
             dtype = sp.float32).reshape((fpa.native_rows, fpa.native_rows))
        self.crf_correction = sp.fromfile(self.crf_correction_file,
             dtype = sp.float32).reshape((fpa.native_columns, fpa.native_columns))
        self.bad = sp.fromfile(self.bad_element_file,
             dtype = sp.int16).reshape((fpa.native_rows, fpa.native_columns))
        self.flat_field = sp.fromfile(self.flat_field_file,
             dtype = sp.float32).reshape((1, fpa.native_rows, fpa.native_columns))
        self.flat_field = self.flat_field[0,:,:]
        _, self.radiometric_calibration, self.radiometric_uncert = \
             sp.loadtxt(self.radiometric_coefficient_file).T

        # Load ghost configuration and construct the matrix
        with open(self.ghost_map_file,'r') as fin:
            ghost_config = json.load(fin)
        self.ghost_matrix = build_ghost_matrix(ghost_config, fpa)
        self.ghost_blur_spectral = ghost_config['blur_spectral']
        self.ghost_blur_spatial = ghost_config['blur_spatial']
             
        basis = envi.open(self.linearity_file+'.hdr').load()
        self.linearity_mu = np.squeeze(basis[0,:])
        self.linearity_mu[np.isnan(self.linearity_mu)] = 0
        self.linearity_evec = np.squeeze(basis[1:,:].T)
        self.linearity_evec[np.isnan(self.linearity_evec)] = 0
        self.linearity_coeffs = envi.open(self.linearity_map_file+'.hdr').load()

@ray.remote
def calibrate_raw(frame, fpa, config):

    # Detector corrections
    frame = subtract_dark(frame, config.dark)
    frame = fix_pedestal(frame, fpa)
    frame = fix_linearity(frame, config.linearity_mu, 
        config.linearity_evec, config.linearity_coeffs)
    frame = frame * config.flat_field
    frame = fix_bad(frame, config.bad, fpa)

    # Optical corrections
    frame = fix_scatter(frame, config.srf_correction, config.crf_correction)
    frame = fix_ghost(frame, fpa, config.ghost_matrix, 
         blur_spatial = config.ghost_blur_spatial, 
         blur_spectral = config.ghost_blur_spectral)

    # Absolute radiometry
    frame = (frame.T * config.radiometric_calibration).T
   
    # Fix OSF
    frame = fix_osf(frame, fpa)

    # Catch NaNs
    frame[sp.logical_not(sp.isfinite(frame))]=0

    # Clip the channels to the appropriate size, if needed
    if config.extract_subframe:
        frame = frame[:,fpa.first_distributed_column:(fpa.last_distributed_column + 1)]
        frame = frame[fpa.first_distributed_row:(fpa.last_distributed_row + 1),:]
        frame = sp.flip(frame, axis=0)

    return frame
   

def main():

    description = "Spectroradiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    default_config = my_directory + '/config/tvac2_config.json'
    parser.add_argument('--config_file', default = default_config)
    parser.add_argument('--dark_file', default = None)
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--maxjobs', type=int, default=30)
    parser.add_argument('input_file', default='')
    parser.add_argument('output_file', default='')
    args = parser.parse_args()

    fpa = FPA(args.config_file)
    config = Config(fpa, args.config_file, args.dark_file)
    ray.init()

    # Set up logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', 
            level=args.level, filename=args.log_file)

    logging.info('Starting calibration')
    lines = 0
    raw = 'Start'

    infile = envi.open(find_header(args.input_file))

    if int(infile.metadata['data type']) == 2:
        dtype = np.int16
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
    lines_analyzed = 0

    with open(args.input_file,'rb') as fin:
        with open(args.output_file,'wb') as fout:

            raw = sp.fromfile(fin, count=nframe, dtype=dtype)
            jobs = []
               
            while len(raw)>0:

                if dtype == np.int16:
                   # left shift by 2 binary digits, 
                   # returning to the 16 bit range.
                   raw = raw * 4

                # Read a frame of data
                if lines_analyzed%10==0:
                    logging.info('Calibrating line '+str(lines_analyzed))
                
                raw = np.array(raw, dtype=sp.float32)

                # All operations take place assuming 480 rows.
                # EMIT avionics only downlink a subset of this data
                # We embed the raw data in a larger frame for analysis.
                frame = raw.reshape((rows,columns))
                if raw.shape[0] < fpa.native_rows:
                    frame = frame_embed(frame, fpa)               

                jobs.append(calibrate_raw.remote(frame, fpa, config))
                lines_analyzed = lines_analyzed + 1

                if len(jobs) == args.maxjobs:
                    
                    # Write to file
                    result = ray.get(jobs)
                    for frame in result:
                        sp.asarray(frame, dtype=sp.float32).tofile(fout)
                    jobs = []
            
                # Read next chunk
                raw = sp.fromfile(fin, count=nframe, dtype=dtype)

            # Do any final jobs
            result = ray.get(jobs)
            for frame in result:
                sp.asarray(frame, dtype=sp.float32).tofile(fout)

    # Form output metadata strings
    wl = config.wl_full.copy()
    fwhm = config.fwhm_full.copy()

    if config.extract_subframe:
        ncolumns = fpa.last_distributed_column - fpa.first_distributed_column + 1
        nchannels = fpa.last_distributed_row - fpa.first_distributed_row + 1
        clip_rows = np.arange(fpa.last_distributed_row, fpa.first_distributed_row-1,-1,dtype=int)
        wl = wl[clip_rows]
        fwhm = fwhm[clip_rows]
    else:
        nchannels, ncolumns = fpa.native_rows, fpa.native_columns

    band_names_string = ','.join(['channel_'+str(i) \
       for i in range(len(wl))])
    fwhm_string =  ','.join([str(w) for w in fwhm])
    wavelength_string = ','.join([str(w) for w in wl])
    
    params = {'lines': lines}
    params.update(**locals())
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_template.format(**params))

    logging.info('Done')


if __name__ == '__main__':

    main()
