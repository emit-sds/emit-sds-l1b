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
#import ray
import pylab as plt

import time

# Import some EMIT-specific functions
my_directory, my_executable = os.path.split(os.path.abspath(__file__))
sys.path.append(my_directory + '/utils/')

from fpa import FPA, frame_embed, frame_extract
from fixbad import fix_bad
from fixosf import fix_osf
from fixlinearity import fix_linearity
from fixscatter import fix_scatter
from fixghost import fix_ghost
from fixelectronicghost import fix_electronic_ghost
from fixghostraster import build_ghost_matrix
from fixghostraster import build_ghost_blur
from pedestal import fix_pedestal
from darksubtract import subtract_dark
from leftshift import left_shift_twice
from emit2dark import bad_flag, dark_from_file
from angread import read_frames, read_frames_metadata


header_template = """ENVI
description = {{AVIRIS-NG calibrated spectral radiance (units: uW nm-1 cm-2 sr-1)}}
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
ang pge input files = {{{input_files_string}}}
ang pge run command = {{{run_command_string}}}
"""

 
def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

     
class Config:

    def __init__(self, fpa, mode):

        # Load calibration file data
        current_mode   = fpa.modes[mode]
        self.radiometric_coefficient_file = current_mode['radiometric_coefficient_file']
        self.flat_field_file = current_mode['flat_field_file']
        self.linearity_file = current_mode['linearity_file']
        self.linearity_map_file = current_mode['linearity_map_file']

        _, self.wl_full, self.fwhm_full = \
             sp.loadtxt(fpa.spectral_calibration_file).T * 1000
        self.srf_correction = sp.fromfile(fpa.srf_correction_file,
             dtype = sp.float32).reshape((fpa.native_rows, fpa.native_rows))
        self.crf_correction = sp.fromfile(fpa.crf_correction_file,
             dtype = sp.float32).reshape((fpa.native_columns, fpa.native_columns))
        self.bad = sp.fromfile(fpa.bad_element_file,
             dtype = sp.int16).reshape((fpa.native_rows, fpa.native_columns))

        #TODO - confirm that 64-bit float is correct
        self.flat_field = sp.fromfile(self.flat_field_file,
             dtype = sp.float32).reshape((2, fpa.native_rows, fpa.native_columns))

        # Take only first band of flat field - second band is the uncertainty
        self.flat_field = self.flat_field[0,:,:]
        self.flat_field[np.logical_not(np.isfinite(self.flat_field))] = 0
        self.radiometric_calibration, self.radiometric_uncert, _ = \
             sp.loadtxt(self.radiometric_coefficient_file).T

        #Truncate to match C code read-in
        self.radiometric_calibration = (1e-7)*np.floor(self.radiometric_calibration*(1e7)) 

        basis = envi.open(self.linearity_file+'.hdr').load()
        self.linearity_mu = np.squeeze(basis[0,:]).copy()
        self.linearity_mu[np.isnan(self.linearity_mu)] = 0
        self.linearity_evec = np.squeeze(basis[1:,:].T).copy()
        self.linearity_evec[np.isnan(self.linearity_evec)] = 0
        self.linearity_coeffs = envi.open(self.linearity_map_file+'.hdr').load()

#@ray.remote
def calibrate_raw_remote(frames, fpa, config):
    return calibrate_raw(frames, fpa, config)

def calibrate_raw(frames, fpa, config):

  if len(frames.shape) == 2:
    frames = np.reshape(frames,(1,frames.shape[0],frames.shape[1]))

  noises = []
  output_frames = []
  for _f in range(frames.shape[0]):
    frame = frames[_f,...]
    noise = -9999

    ## Don't calibrate a bad frame
    if not np.all(frame <= bad_flag):

        # Left shift, returning to the 16 bit range.
        if hasattr(fpa,'left_shift_twice') and fpa.left_shift_twice:
           frame = left_shift_twice(frame)

        # Dark state subtraction
        frame = subtract_dark(frame, config.dark)
       
        ## Delete telemetry
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

        # Electronic ghost
        frame = fix_electronic_ghost(frame, fpa.eghost_samples_per_panel, np.array(fpa.eghost_template),
                                     fpa.eghost_panel_correction, fpa.eghost_panel_multipliers)
     
        ##frame = fix_linearity(frame, config.linearity_mu, 
        ##    config.linearity_evec, config.linearity_coeffs)
        # frame[frame<0] = 0
        # frame[frame>=65536] = 65535

        temp_flatfield = config.flat_field.copy()
        temp_flatfield[temp_flatfield==0] = 1.0
        temp_flatfield[temp_flatfield>2.0] = 1.0
        temp_flatfield[temp_flatfield<0.5] = 1.0

        frame = frame * temp_flatfield
        
        #frame = frame * config.flat_field
        
        # Fix bad pixels, and any nonfinite results from the previous
        # operations
        flagged = np.logical_not(np.isfinite(frame))
        frame[flagged] = 0
        bad = config.bad.copy()
        bad[flagged] = -1
        frame = fix_bad(frame, bad, fpa)

        frame[frame<-50] = 0.0
        frame[frame>1e9] = 0.0
        frame[np.logical_not(np.isfinite(frame))] = 0.0
        
        # Optical corrections
        frame = fix_scatter(frame, config.srf_correction, config.crf_correction)
        
        # Absolute radiometry
        frame = (frame.T * config.radiometric_calibration).T
       
        # Fix OSF
        #frame = fix_osf(frame, fpa)
        
        # Catch NaNs
        frame[sp.logical_not(sp.isfinite(frame))]=0

    # Clip the channels to the appropriate size, if needed
    if fpa.extract_subframe:
        frame = frame[:,fpa.first_distributed_column:(fpa.last_distributed_column + 1)]
        frame = frame[fpa.first_distributed_row:(fpa.last_distributed_row + 1),:]
        frame = sp.flip(frame, axis=0)


    output_frames.append(frame)
    noises.append(noise)

  # Replace all bad data flags with -9999
  output_frames = np.stack(output_frames)
  output_frames[output_frames<=(bad_flag+1e-6)] = np.nan

  noises = np.array(noises)
  if np.sum(noises != -9999) > 0:
    noises = np.nanmedian(noises[noises != -9999])
  else:
    noises = -9999

  # Co-add
  output_frames = np.nanmean(output_frames,axis=0)
  output_frames[np.isnan(output_frames)] = -9999

  return output_frames, noises
   

def main():

    description = "Spectroradiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--mode', default = 'default')
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--max_jobs', type=int, default=40)
    parser.add_argument('--debug_mode', action='store_true')
    parser.add_argument('input_file', default='')
    parser.add_argument('config_file', default='')
    parser.add_argument('output_file', default='')
    args = parser.parse_args()

    start_time = time.time()
    
    fpa = FPA(args.config_file)
    config = Config(fpa, args.mode)
    
    #Find binfac
    binfac_file = args.input_file + '.binfac'
    if os.path.isfile(binfac_file) is False:
        logging.error(f'binfac file not found at expected location: {binfac_file}')
        raise ValueError('Binfac file not found - see log for details')
    binfac = int(np.genfromtxt(binfac_file))
    
    # Set up logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', 
            level=args.level, filename=args.log_file)
    
    logging.info('Starting calibration')
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
    
    rows = int(infile.metadata['bands']) - 1 # extra band is metadata
    columns = int(infile.metadata['samples'])
    lines_analyzed = 0
    noises = []
    
    # Read metadata from RAW ang file
    logging.debug('Reading metadata')
    frame_meta, num_read, frame_obcv = read_frames_metadata(args.input_file, 500000, rows, columns, 0)
    
    dark_frame_idxs = np.where(frame_obcv == 2)[0]
    science_frame_idxs = np.where(frame_obcv[dark_frame_idxs[-1]+1:])[0] + dark_frame_idxs[-1] + 1
    
    #*** need to read this in fro GLT for consistency with c code
    
    dark_frame_start_idx = dark_frame_idxs[fpa.dark_margin] # Trim to make sure the shutter transition isn't in the dark
    num_dark_frames = 600 #Matches C code directly
    #dark_frame_idxs[-1*fpa.dark_margin]-dark_frame_start_idx
    # ASSUMES THAT OBCV only gives darks at the beginning
    
    logging.debug(f'Found {len(dark_frame_idxs)} dark frames and {len(science_frame_idxs)} science frames')
    
    if np.all(science_frame_idxs - science_frame_idxs[0] == np.arange(len(science_frame_idxs))) is False:
        logging.error('Science frames are not contiguous, cannot proceed')
        raise AttributeError('Science frames are not contiguous')

    # Read dark
    dark_frames, _, _, _ = read_frames(args.input_file, 599, fpa.native_rows, fpa.native_columns, 101)
    config.dark = np.mean(dark_frames,axis=0) #np.median(dark_frames,axis=0)
    config.dark_std = np.std(dark_frames,axis=0)
    del dark_frames
    logging.debug('Dark read complete, beginning calibration')
    #ray.init()
    #fpa_id = ray.put(fpa)

    jobs = []
    if args.debug_mode:
        result = []
    for sc_idx in range(science_frame_idxs[0], science_frame_idxs[0] + len(science_frame_idxs), binfac):
        if sc_idx + binfac >= science_frame_idxs[-1]:
            break
    
        frames, frame_meta, num_read, frame_obcv = read_frames(args.input_file, binfac, fpa.native_rows, fpa.native_columns, sc_idx)
    
        if lines_analyzed%10==0:
            logging.info('Calibrating line '+str(lines_analyzed))
    
        if args.debug_mode:
            result.append(calibrate_raw(frames, fpa, config))
        else:
            jobs.append(calibrate_raw_remote.remote(frames, fpa_id, config))
        lines_analyzed = lines_analyzed + 1
    
    num_output_lines = 0
    with open(args.output_file,'wb') as fout:
        # Do any final jobs
        if args.debug_mode is False:
            result = ray.get(jobs)
        for frame, noise in result:
            sp.asarray(frame, dtype=sp.float32).tofile(fout)
            noises.append(noise)
            num_output_lines += 1
    
    # Form output metadata strings
    wl = config.wl_full.copy()
    fwhm = config.fwhm_full.copy()
    
    if fpa.extract_subframe:
        ncolumns = fpa.last_distributed_column - fpa.first_distributed_column + 1
        nchannels = fpa.last_distributed_row - fpa.first_distributed_row + 1
        clip_rows = np.arange(fpa.last_distributed_row, fpa.first_distributed_row-1,-1,dtype=int)
        wl = wl[::-1][clip_rows]
        fwhm = fwhm[::-1][clip_rows]
    else:
        nchannels, ncolumns = fpa.native_rows, fpa.native_columns
    
    band_names_string = ','.join(['channel_'+str(i) \
       for i in range(len(wl))])
    fwhm_string =  ','.join([str(w) for w in fwhm])
    wavelength_string = ','.join([str(w) for w in wl])
    
    params = {}
    params['masked_pixel_noise'] = np.nanmedian(np.array(noises))
    params['run_command_string'] = ' '.join(sys.argv)
    params['input_files_string'] = ''
    for var in dir(fpa):
       if var.endswith('_file'):
          params['input_files_string'] = params['input_files_string'] + \
             ' %s=%s'%(var,getattr(fpa,var))
    params['lines'] =  num_output_lines
    
    params.update(**locals())
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_template.format(**params))
    
    end_time = time.time()
    logging.info(f'Completed {len(science_frame_idxs)} frames in {end_time-start_time} seconds')
    logging.info('Done')


if __name__ == '__main__':

    main()

