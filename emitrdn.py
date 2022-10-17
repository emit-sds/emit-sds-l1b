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
import pylab as plt

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
from fixghostraster import build_ghost_blur
from pedestal import fix_pedestal
from darksubtract import subtract_dark
from leftshift import left_shift_twice
from emit2dark import bad_flag, dark_from_file


header_template = """ENVI
description = {{EMIT L1B calibrated spectral radiance (units: uW nm-1 cm-2 sr-1)}}
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
band names = {{{band_names_string}}}
masked pixel noise = {masked_pixel_noise}
emit pge input files = {{{input_files_string}}}
emit pge run command = {{{run_command_string}}}
flip horizontal  = {flip_horizontal}
"""


replaced_header_template = """ENVI
description = {{EMIT replaced channels}}
samples = {ncolumns}
lines = {lines}
bands = {nreplacedchannels}
header offset = 0
file type = ENVI Standard
data type = 1
interleave = bil
byte order = 0
flip horizontal  = {flip_horizontal}
"""


 
def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

     
class Config:

    def __init__(self, fpa, dark_file, mode):

        # Load calibration file data
        current_mode   = fpa.modes[mode]
        self.radiometric_coefficient_file = current_mode['radiometric_coefficient_file']
        self.flat_field_file = current_mode['flat_field_file']
        self.linearity_file = current_mode['linearity_file']
        self.linearity_map_file = current_mode['linearity_map_file']

        self.dark_frame_file = dark_file
        self.dark, self.dark_std = dark_from_file(self.dark_frame_file)
  
        # Move this outside, to the main function
        if hasattr(fpa,'left_shift_twice') and fpa.left_shift_twice:
           # left shift, returning to the 16 bit range.
           self.dark = left_shift_twice(self.dark)

        _, self.wl_full, self.fwhm_full = \
             sp.loadtxt(fpa.spectral_calibration_file).T * 1000
        self.srf_correction = sp.fromfile(fpa.srf_correction_file,
             dtype = sp.float32).reshape((fpa.native_rows, fpa.native_rows))
        self.crf_correction = sp.fromfile(fpa.crf_correction_file,
             dtype = sp.float32).reshape((fpa.native_columns, fpa.native_columns))
        self.bad = sp.fromfile(fpa.bad_element_file,
             dtype = sp.int16).reshape((fpa.native_rows, fpa.native_columns))
        self.flat_field = sp.fromfile(self.flat_field_file,
             dtype = sp.float32).reshape((1, fpa.native_rows, fpa.native_columns))
        self.flat_field = self.flat_field[0,:,:]
        self.flat_field[np.logical_not(np.isfinite(self.flat_field))] = 0
        _, self.radiometric_calibration, self.radiometric_uncert = \
             sp.loadtxt(self.radiometric_coefficient_file).T

        # zero offset perturbation
        self.zero_offset = np.zeros((fpa.native_rows, fpa.native_columns))
        if hasattr(fpa, 'zero_offset_file'):
            self.zero_offset = sp.fromfile(fpa.zero_offset_file,
                dtype=sp.float32).reshape((1, fpa.native_rows, fpa.native_columns))


        # Load ghost configuration and construct the matrix
        with open(fpa.ghost_map_file,'r') as fin:
            ghost_config = json.load(fin)
        self.ghost_matrix = build_ghost_matrix(ghost_config, fpa)
        self.ghost_blur = build_ghost_blur(ghost_config, fpa)
        self.ghost_center = ghost_config['center']
             
        basis = envi.open(self.linearity_file+'.hdr').load()
        self.linearity_mu = np.squeeze(basis[0,:])
        self.linearity_mu[np.isnan(self.linearity_mu)] = 0
        self.linearity_evec = np.squeeze(basis[1:,:].T)
        self.linearity_evec[np.isnan(self.linearity_evec)] = 0
        self.linearity_coeffs = envi.open(self.linearity_map_file+'.hdr').load()

@ray.remote
def calibrate_raw(frame, fpa, config):

    saturated = np.ones(frame.shape)<0 # False
    bad = config.bad.copy()

    # Don't calibrate a bad frame
    if not np.all(frame <= bad_flag):

        # Left shift, returning to the 16 bit range.
        if hasattr(fpa,'left_shift_twice') and fpa.left_shift_twice:
           frame = left_shift_twice(frame)

        # Test for saturation
        if hasattr(fpa,'saturation_DN'):
            saturated = frame>fpa.saturation_DN
 
        # Dark state subtraction
        frame = subtract_dark(frame, config.dark)
        frame = frame - config.zero_offset
       
        # Delete telemetry
        if hasattr(fpa,'ignore_first_row') and fpa.ignore_first_row:
           frame[0,:] = frame[1,:]
        
        # Raw noise calculation
        if hasattr(fpa,'masked_columns'):
            noise = np.nanmedian(np.std(frame[:,fpa.masked_columns],axis=0))
        elif hasattr(fpa,'masked_rows'):
            noise = np.nanmedian(np.std(frame[fpa.masked_rows,:],axis=1))
        else:
            noise = -1 

        # Detector corrections
        frame = fix_pedestal(frame, fpa)
     
        frame = fix_linearity(frame, config.linearity_mu, 
            config.linearity_evec, config.linearity_coeffs)
        frame = frame * config.flat_field
        
        # Fix bad pixels, saturated pixels, and any nonfinite 
        # results from the previous operations
        flagged = np.logical_or(saturated, np.logical_not(np.isfinite(frame)))
        frame[flagged] = 0
        bad[flagged] = -1
        frame = fix_bad(frame, bad, fpa)
        
        # Optical corrections
        frame = fix_scatter(frame, config.srf_correction, config.crf_correction)
        frame = fix_ghost(frame, fpa, config.ghost_matrix, 
             blur = config.ghost_blur, center = config.ghost_center)
        
        # Absolute radiometry
        frame = (frame.T * config.radiometric_calibration).T
       
        # Fix OSF
        frame = fix_osf(frame, fpa)
        
        # Catch NaNs
        frame[sp.logical_not(sp.isfinite(frame))]=0

    else:
        noise = -9999

    if fpa.extract_subframe:

        # Clip the radiance data to the appropriate size
        frame = frame[:,fpa.first_distributed_column:(fpa.last_distributed_column + 1)]
        frame = frame[fpa.first_distributed_row:(fpa.last_distributed_row + 1),:]
        frame = sp.flip(frame, axis=0)

        # Clip the replaced channel mask
        bad = bad[:,fpa.first_distributed_column:(fpa.last_distributed_column + 1)]
        bad = bad[fpa.first_distributed_row:(fpa.last_distributed_row + 1),:]
        bad = sp.flip(bad, axis=0)

    # Mirror image
    if hasattr(fpa, 'flip_horizontal') and fpa.flip_horizontal:
        bad = sp.flip(bad, axis=1)
        frame = sp.flip(frame, axis=1)

    # Replace all bad data flags with -9999
    cleanframe = frame.copy()
    cleanframe[frame<=(bad_flag+1e-6)] = -9999

    return cleanframe, noise, np.packbits(bad, axis=0) 
   

def main():

    description = "Spectroradiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--mode', default = 'default')
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--max_jobs', type=int, default=40)
    parser.add_argument('input_file', default='')
    parser.add_argument('dark_file', default = None)
    parser.add_argument('config_file', default='')
    parser.add_argument('output_file', default='')
    parser.add_argument('output_replaced', default='')
    args = parser.parse_args()

    fpa = FPA(args.config_file)
    config = Config(fpa, args.dark_file, args.mode)
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
    lines_analyzed = 0
    noises = []

    with open(args.input_file,'rb') as fin:
       with open(args.output_file,'wb') as fout:
          with open(args.output_replaced,'wb') as foutreplace:

            raw = sp.fromfile(fin, count=nframe, dtype=dtype)
            jobs = []
               
            while len(raw)>0:

                # Read a frame of data 
                raw = np.array(raw, dtype=sp.float32)
                frame = raw.reshape((rows,columns))

                if lines_analyzed%10==0:
                    logging.info('Calibrating line '+str(lines_analyzed))

                jobs.append(calibrate_raw.remote(frame, fpa, config))
                lines_analyzed = lines_analyzed + 1

                if len(jobs) == args.max_jobs:
                    
                    # Write to file
                    result = ray.get(jobs)
                    for frame, noise, bad in result:
                        np.asarray(frame, dtype=sp.float32).tofile(fout)
                        np.asarray(bad, dtype=sp.uint8).tofile(foutreplace)
                        noises.append(noise)
                    jobs = []
            
                # Read next chunk
                raw = sp.fromfile(fin, count=nframe, dtype=dtype)

            # Do any final jobs
            result = ray.get(jobs)
            for frame, noise, bad in result:
                sp.asarray(frame, dtype=sp.float32).tofile(fout)
                np.asarray(bad, dtype=sp.uint8).tofile(foutreplace)
                noises.append(noise)

    # Form output metadata strings
    wl = config.wl_full.copy()
    fwhm = config.fwhm_full.copy()

    if fpa.extract_subframe:
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
    
    # Place all calibration parameters in header metadata
    params = {'lines': lines}
    params['masked_pixel_noise'] = np.nanmedian(np.array(noises))
    params['run_command_string'] = ' '.join(sys.argv)
    params['input_files_string'] = ' dark_file='+args.dark_file
    for var in dir(fpa):
       if var.endswith('_file'):
          params['input_files_string'] = params['input_files_string'] + \
             ' %s=%s'%(var,getattr(fpa,var))

    flip_horizontal = None
    if hasattr(fpa, 'flip_horizontal') and fpa.flip_horizontal:
        flip_horizontal = 1
    else:
        flip_horizontal = 0

    # Write the header
    params.update(**locals())
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_template.format(**params))


    # Output the header file for the replaced pixel image
    nreplacedchannels = bad.shape[0]
    params = {'lines': lines}
    params.update(**locals())
    with open(args.output_replaced+'.hdr','w') as fout:
        fout.write(replaced_header_template.format(**params))
    logging.info('Done')


if __name__ == '__main__':

    main()
