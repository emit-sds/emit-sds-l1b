#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import scipy as sp
import numpy as np
from spectral.io import envi
from datetime import datetime, timezone
from scipy import linalg, polyfit, polyval
import json
import logging
import argparse
import multiprocessing

software_version = "0.0.1"
product_version = "0.0.1"
documentation_version = "EMIT L1B ATBD v1.0"


header_template = """ENVI
description = {{Calibrated Radiance, microWatts per (steradian nanometer [centemeter squared])}}
samples = {columns}
lines = {lines}
bands = {channels}
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
byte order = 0
wavelength units = Nanometers
wavelength = {{{wavelength_string}}}
fwhm = {{{fwhm_string}}}
band names = {{{band_names_string}}}
emit_acquisition start time = YYYYMMDDTHHMMSS
emit_acquisition stop time = YYYYMMDDTHHMMSS
emit_pge_name = emit-sds-l1b/l1a_rdn
emit_pge_version = {software_version}
emit_pge_input files = {{
  dn_file = {input_file} 
  radiometric_coefficient_file = {radiometric_coefficient_file} 
  spectral_calibration_file = {spectral_calibration_file} 
  srf_correction_file = {srf_correction_file} 
  crf_correction_file = {crf_correction_file} 
  bad_element_file = {bad_element_file} 
  flat_field_file = {flat_field_file} 
  radiometric_coefficient_file = {radiometric_coefficient_file} 
  linearity_file = {linearity_file}
  }}
emit software build version = {software_version}
emit documentation version = {documentation_version}
emit data product creation time = {creation_time} 
emit data product version = {product_version}
emit orbit correction performed = false """


      
class Config:

    def __init__(self, filename, input_file='', output_file=''):

        # load calibration file data
        ok = 'OK!'
        with open(filename,'r') as fin:
         self.__dict__ = json.load(fin)
        try:
           self.dark, _ = sp.fromfile(self.dark_frame_file,
                dtype = sp.float32).reshape((2,self.channels, self.columns))
           _, self.wl, self.fwhm = \
                sp.loadtxt(self.spectral_calibration_file).T * 1000
           self.srf_correction = sp.fromfile(self.srf_correction_file,
                dtype = sp.float32).reshape((self.channels, self.channels))
           self.crf_correction = sp.fromfile(self.crf_correction_file,
                dtype = sp.float32).reshape((self.columns, self.columns))
           self.bad = sp.fromfile(self.bad_element_file,
                dtype = sp.uint16).reshape((self.channels, self.columns))
           self.flat_field, _ = sp.fromfile(self.flat_field_file,
                dtype = sp.float32).reshape((2,self.channels, self.columns))
           self.radiometric_calibration, _, _ = \
                sp.loadtxt(self.radiometric_coefficient_file).T
           self.linearity = sp.fromfile(self.linearity_file, 
                dtype=sp.uint16).reshape((65536,))
        except ValueError:
            logging.error('Incorrect file size for calibration data')
        except AttributeError:
            logging.error('One or more missing calibration files')

        # Check for NaNs in calibration data
        for M in [self.dark, self.wl, self.srf_correction, 
                self.crf_correction, self.bad, self.flat_field,
                self.radiometric_calibration, self.linearity]:
            if (np.logical_not(sp.isfinite(M))).sum() > 0:
                logging.error('Non-finite values in calibration data')

        # Truncate flat field values, if needed
        if self.flat_field_limits is not None:
           lo, hi = self.flat_field_limits
           self.flat_field[self.flat_field < lo] = lo
           self.flat_field[self.flat_field > hi] = hi 

        # size of regular frame and raw frame (with header)
        self.frame_shape = (self.channels, self.columns)
        self.nframe = sp.prod(self.frame_shape)
        self.raw_shape = (self.channels + self.header_channels, self.columns)
        self.nraw = sp.prod(self.raw_shape)

        # form output metadata strings
        self.band_names_string = ','.join(['channel_'+str(i) \
                for i in range(len(self.wl))])
        self.fwhm_string =  ','.join([str(w) for w in self.fwhm])
        self.wavelength_string = ','.join([str(w) for w in self.wl])

        # find the input files
        if len(input_file)>0:
            self.input_file = input_file
        if len(output_file)>0:
            self.output_file = output_file

        # identify input file header
        if self.input_file.endswith('.img'):
            self.input_header = self.input_file.replace('.img','.hdr') 
        else:
            self.input_header = self.input_file + '.hdr'

        # identify output file header
        if self.output_file.endswith('.img'):
            self.output_header = self.output_file.replace('.img','.hdr') 
        else:
            self.output_header = self.output_file + '.hdr'

        # read acquisition start and stop from the L1a
        meta = envi.read_envi_header(self.input_header)
        if 'emit acquisition start time' in meta:
            self.emit_acquisition_start_time = \
                meta[emit_acquisition_start_time]
        else:
            self.emit_acquisition_time = 'YYYYMMDDTHHMMSS'
        if 'emit acquisition stop time' in meta:
            self.emit_acquisition_stop_time = \
                meta[emit_acquisition_stop_time]
        else:
            self.emit_acquisition_time = 'YYYYMMDDTHHMMSS'


def correct_pedestal_shift(frame, config):
    mean_dark = frame[config.dark_channels,:].mean(axis=0)
    return frame - mean_dark


def infer_bad(frame, col, config):
    '''Infer the value of a bad pixel'''
    clean = sp.where(np.logical_not(config.bad).all(axis=1))[0]
    bad = config.bad[:,col]
    sa = frame[clean,:].T @ frame[clean, col]
    norms = linalg.norm(frame[clean,:], axis=0).T
    sa = sa / (norms * norms[col])
    sa[col] = 9e99
    best = sp.argmin(sa)
    p = polyfit(frame[clean, best], frame[clean, col],1)
    new = frame[:,col]
    new[bad] = polyval(p, frame[bad, best])
    return new 

    
def fix_bad(frame, config):
    fixed = frame.copy()
    for col in sp.nonzero(config.bad.any(axis=0))[0]:
        fixed[:,col] = infer_bad(frame, col, config)
    return fixed


def subtract_dark(frame, dark):
    return frame - dark


def correct_spatial_resp(frame, crf_correction):
    scratch = sp.zeros(frame.shape)
    for i in range(frame.shape[0]):
        scratch[i,:] = frame[i,:] @ crf_correction 
    return scratch


def correct_spectral_resp(frame, srf_correction):
    scratch = sp.zeros(frame.shape)
    for i in range(frame.shape[1]):
        scratch[:,i] = frame[:,i].T @  srf_correction
    return scratch


def detector_corrections(frame, config):
    return frame


def correct_panel_ghost(frame, config):

    pg_template = sp.array([config.pg_template])
    ntemplate = len(config.pg_template)

    panel1 = sp.arange(config.panel_width)
    panel2 = sp.arange(config.panel_width,(2*config.panel_width))
    panel3 = sp.arange((2*config.panel_width),(3*config.panel_width))
    panel4 = sp.arange((3*config.panel_width),(4*config.panel_width))
 
    avg1 = frame[:,panel1].mean(axis=1)[:,sp.newaxis]
    avg2 = frame[:,panel2].mean(axis=1)[:,sp.newaxis]
    avg3 = frame[:,panel3].mean(axis=1)[:,sp.newaxis]
    avg4 = frame[:,panel4].mean(axis=1)[:,sp.newaxis]
  
    c1 = frame[:,panel1];
    c2 = frame[:,panel2];
    c3 = frame[:,panel3];
    c4 = frame[:,panel4];       
 
    coef1 = config.panel_ghost_correction*(c2+c3+c4);
    coef2 = config.panel_ghost_correction*(c1+c3+c4);
    coef3 = config.panel_ghost_correction*(c1+c2+c4);
    coef4 = config.panel_ghost_correction*(c1+c2+c3);       

    coef1[:,:ntemplate] = 1.6 * (avg2+avg3+avg4) @ pg_template
    coef2[:,:ntemplate] = 1.6 * (avg1+avg3+avg4) @ pg_template
    coef3[:,:ntemplate] = (avg1+avg2+avg4)@ pg_template
    coef4[:,:ntemplate] = (avg1+avg2+avg3)@ pg_template
            
    new = sp.zeros(frame.shape)
    new[:,panel1] = frame[:,panel1] + coef1;
    new[:,panel2] = frame[:,panel2] + coef2;
    new[:,panel3] = frame[:,panel3] + coef3;
    new[:,panel4] = frame[:,panel4] + coef4;

    return new


def main():

    description = "Radiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config_file')
    parser.add_argument('input_file', nargs='?', default='')
    parser.add_argument('output_file', nargs='?', default='')
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s {version}'.format(version=software_version))
    args = parser.parse_args()

    config = Config(args.config_file, args.input_file, args.output_file)
        
    logging.basicConfig(format='%(message)s', level=args.level)
    
    
    logging.info('Starting calibration')
    lines = 0
    raw = 'Start'

    print(config.input_file)
    with open(config.input_file,'rb') as fin:
        with open(config.output_file,'wb') as fout:

            raw = sp.fromfile(fin, count=config.nraw, dtype=sp.uint16)
            while len(raw)>0:

                # Read a frame of data
                if lines%10==0:
                    logging.info('Calibrating line '+str(lines))
                
                raw = raw.reshape(config.raw_shape)
                header = raw[:config.header_channels, :]
                frame  = raw[config.header_channels:,:]
                
                # Detector corrections
                frame = subtract_dark(frame, config.dark)
                frame = correct_pedestal_shift(frame, config)
                frame = correct_panel_ghost(frame, config) 
                frame = (frame.T * config.radiometric_calibration).T
                frame = frame * config.flat_field
                frame = fix_bad(frame, config)
                frame = correct_spectral_resp(frame, config.srf_correction)
                frame = correct_spatial_resp(frame, config.crf_correction)
   
                # Reverse channels, catch NaNs, and write
                frame[sp.logical_not(sp.isfinite(frame))]=0
                if config.reverse_channels:
                    frame = sp.flip(frame, axis=0)
                sp.asarray(frame, dtype=sp.float32).tofile(fout)
                lines = lines + 1
            
                # read next chunk
                raw = sp.fromfile(fin, count=config.nraw, dtype=sp.uint16)

    params = {'lines': lines,
     'creation_time': datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}
    params.update(globals())
    params.update(config.__dict__)
    with open(config.output_header,'w') as fout:
        fout.write(header_template.format(**params))

    logging.info('Done')

if __name__ == '__main__':

    main()
