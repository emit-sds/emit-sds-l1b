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

# Import some EMIT-specific functions
my_directory = os.path.abspath(__file__)
sys.path.append(my_directory + 'utlis/')

from emit_fpa import native_rows, frame_embed, frame_extract
from fixbad import fix_bad 
from fixlinearity import fix_frame
from fixscatter import fix_scatter
from fixghost import fix_ghost
from pedestal import fix_pedestal

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
band names = {{{band_names_string}}}"""

      
class Config:

    def __init__(self, filename, input_file='', output_file='', dark_file=None):

        # Load calibration file data
        with open(filename,'r') as fin:
         self.__dict__ = json.load(fin)
        try:
           if dark_file is not None:
               self.dark_frame_file = dark_file
           self.dark = darkfromfile(self.dark_frame_file)
           _, self.wl_full, self.fwhm_full = \
                sp.loadtxt(self.spectral_calibration_file).T * 1000
           self.srf_correction = sp.fromfile(self.srf_correction_file,
                dtype = sp.float32).reshape((self.channels_raw, 
                    self.channels_raw))
           self.crf_correction = sp.fromfile(self.crf_correction_file,
                dtype = sp.float32).reshape((self.columns_raw, 
                    self.columns_raw))
           self.bad = sp.fromfile(self.bad_element_file,
                dtype = sp.uint16).reshape((self.channels_raw, 
                    self.columns_raw))
           self.flat_field = sp.fromfile(self.flat_field_file,
                dtype = sp.float32).reshape((2, self.channels_raw, 
                    self.columns_raw))[0,:,:]
           self.radiometric_calibration, _, _ = \
                sp.loadtxt(self.radiometric_coefficient_file).T
           self.linearity = sp.fromfile(self.linearity_file, 
                dtype=sp.uint16).reshape((65536,))
        except ValueError:
            logging.error('Incorrect file size for calibration data')
        except AttributeError:
            logging.error('One or more missing calibration files')

        if self.channels != (self.channels_raw - \
            (self.channels_masked[0] + self.channels_masked[1])):
            raise ValueError('channel mask inconsistent with total channels')

        if self.columns != (self.columns_raw - \
            (self.columns_masked[0] + self.columns_masked[1])):
            raise ValueError('column mask inconsistent with total channels')

        self.wl = np.array([w for w in \
            self.wl_full[self.channels_masked[0]:-self.channels_masked[1]]])
        self.fwhm = np.array([w for w in \
            self.fwhm_full[self.channels_masked[0]:-self.channels_masked[1]]])

        # Check for NaNs in calibration data
        for name in ['dark', 'wl_full', 'srf_correction', 
                'crf_correction', 'bad', 'flat_field',
                'radiometric_calibration','linearity']:
            obj = getattr(self, name)
            invalid  = np.logical_not(sp.isfinite(obj))
            if invalid.sum() > 0:
                msg='Replacing %i non-finite values in %s' 
                logging.warning(msg % (invalid.sum(),name))
            obj[invalid]=0 

        # Size of regular frame and raw frame (with header)
        self.frame_shape = (self.channels, self.columns)
        self.nframe = sp.prod(self.frame_shape)
        self.raw_shape = (self.channels_raw + self.header_channels, self.columns_raw)
        self.nraw = sp.prod(self.raw_shape)

        # Form output metadata strings
        self.band_names_string = ','.join(['channel_'+str(i) \
                for i in range(len(self.wl))])
        self.fwhm_string =  ','.join([str(w) for w in self.fwhm])
        self.wavelength_string = ','.join([str(w) for w in self.wl])

        # Clean channels have no bad elements
        self.clean = sp.where(np.logical_not(self.bad).all(axis=1))[0]
        logging.warning(str(len(self.clean))+' clean channels')

        # Find the input files
        if len(input_file)>0:
            self.input_file = input_file
        if len(output_file)>0:
            self.output_file = output_file

        # Identify input file header
        if self.input_file.endswith('.img'):
            self.input_header = self.input_file.replace('.img','.hdr') 
        else:
            self.input_header = self.input_file + '.hdr'

        # Identify output file header
        if self.output_file.endswith('.img'):
            self.output_header = self.output_file.replace('.img','.hdr') 
        else:
            self.output_header = self.output_file + '.hdr'




def main():

    description = "Spectroradiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    my_directory, my_file = os.path.split(os.path.realpath(__file___))
    default_config = my_directory + '/configs/tvac2_config.json'
    parser.add_argument('--config_file', default = default_config)
    parser.add_argument('--dark_file', default = None)
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('input_file', nargs='?', default='')
    parser.add_argument('output_file', nargs='?', default='')
    args = parser.parse_args()

    config = Config(args.config_file, args.input_file, 
        args.output_file, args.dark_file)

    # Set up logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.level, filename=args.log_file)

    logging.info('Starting calibration')
    lines = 0
    raw = 'Start'

    with open(config.input_file,'rb') as fin:
        with open(config.output_file,'wb') as fout:

            raw = sp.fromfile(fin, count=config.nraw, dtype=sp.int16)
            while len(raw)>0:

                # Read a frame of data
                if lines%10==0:
                    logging.info('Calibrating line '+str(lines))
                
                raw = np.array(raw, dtype=sp.float32)
                raw = raw.reshape(config.raw_shape)
                header = raw[:config.header_channels, :]
                frame  = raw[config.header_channels:, :]
                if frame.shape[0] < native_rows:
                    frame = frame_embed(frame)
                
                # Detector corrections
                frame = subtract_dark(frame, config.dark)
                frame = correct_pedestal_shift(frame)
                frame = fix_bad(frame, config.badmap)
                frame = frame * config.flat_field

                # Optical corrections
                frame = fix_scatter(frame, config.spectral_correction, config.spatial_correction)
                frame = fix_ghost(frame, config.ghostmap)

                # Absolute radiometry
                frame = (frame.T * config.radiometric_calibration).T
   
                # Reverse channels, catch NaNs, and write
                frame[sp.logical_not(sp.isfinite(frame))]=0
                if config.reverse_channels:
                    frame = sp.flip(frame, axis=0)

                # Clip the channels to the appropriate size
                if config.extract_subframe:
                    frame = frame_extract(frame, columns = True)
                sp.asarray(frame, dtype=sp.float32).tofile(fout)
                lines = lines + 1
            
                # Read next chunk
                raw = sp.fromfile(fin, count=config.nraw, dtype=sp.int16)

    params = {'lines': lines}
    params.update(globals())
    params.update(config.__dict__)
    with open(config.output_header,'w') as fout:
        fout.write(header_template.format(**params))

    logging.info('Done')


if __name__ == '__main__':

    main()
