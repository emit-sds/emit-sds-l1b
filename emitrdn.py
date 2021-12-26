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
my_directory, my_executable = os.path.split(os.path.abspath(__file__))
sys.path.append(my_directory + '/utils/')
print(sys.path)

from emit_fpa import native_rows, native_columns, frame_embed, frame_extract
from fixbad import fix_bad 
from fixlinearity import fix_frame
from fixscatter import fix_scatter
from fixghost import fix_ghost
from pedestal import fix_pedestal
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

    def __init__(self, filename, input_file='', output_file='', dark_file=None):

        # Load calibration file data
        with open(filename,'r') as fin:
         self.__dict__ = json.load(fin)
        
        # Adjust local filepaths where needed
        for fi in ['spectral_calibration_file','srf_correction_file',
                   'crf_correction_file','linearity_file',
                   'radiometric_coefficient_file',
                   'bad_element_file','flat_field_file']:
            path = getattr(self,fi)
            if path[0] != '/':
                path = os.path.join(my_directory, path)
                setattr(self,fi,path)

        try:
           if dark_file is not None:
               self.dark_frame_file = dark_file
           self.dark, self.dark_std = dark_from_file(self.dark_frame_file)
           _, self.wl_full, self.fwhm_full = \
                sp.loadtxt(self.spectral_calibration_file).T * 1000
           self.srf_correction = sp.fromfile(self.srf_correction_file,
                dtype = sp.float32).reshape((native_rows, native_rows))
           self.crf_correction = sp.fromfile(self.crf_correction_file,
                dtype = sp.float32).reshape((native_columns, native_columns))
           self.bad = sp.fromfile(self.bad_element_file,
                dtype = sp.uint16).reshape((native_rows, native_columns))
           self.flat_field = sp.fromfile(self.flat_field_file,
                dtype = sp.float32).reshape((2, native_rows, native_columns))
           self.flat_field = self.flat_field[0,:,:]
           _, self.radiometric_calibration, self.radiometric_uncert = \
                sp.loadtxt(self.radiometric_coefficient_file).T
           self.linearity = sp.fromfile(self.linearity_file, 
                dtype=sp.uint16, count=65536).reshape((65536,))
        except ValueError:
            logging.error('Incorrect file size for calibration data')
        except AttributeError:
            logging.error('One or more missing calibration files')

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
    default_config = my_directory + '/configs/tvac2_config.json'
    parser.add_argument('--config_file', default = default_config)
    parser.add_argument('--dark_file', default = None)
    parser.add_argument('--level', default='DEBUG',
            help='verbosity level: INFO, ERROR, or DEBUG')
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('input_file', default='')
    parser.add_argument('output_file', default='')
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

    with open(config.input_file,'rb') as fin:
        with open(config.output_file,'wb') as fout:

            raw = sp.fromfile(fin, count=nframe, dtype=dtype)

            if dtype == np.int16:
               # left shift by 2 binary digits, 
               # returning to the 16 bit range.
               raw = raw * 4
               
            while len(raw)>0:

                # Read a frame of data
                if lines%10==0:
                    logging.info('Calibrating line '+str(lines))
                
                raw = np.array(raw, dtype=sp.float32)
                frame = raw.reshape((rows,columns))

                # All operations take place assuming 480 rows.
                # EMIT avionics only downlink a subset of this data
                if raw.shape[0] < native_rows:
                    frame = frame_embed(frame)
                
                # Detector corrections
                frame = subtract_dark(frame, config.dark)
                frame = correct_pedestal_shift(frame)
                frame = fix_bad(frame, config.badmap)
                frame = frame * config.flat_field

                # Optical corrections
                frame = fix_scatter(frame, 
                    config.spectral_correction, 
                    config.spatial_correction)
                frame = fix_ghost(frame, config.ghostmap)

                # Absolute radiometry
                frame = (frame.T * config.radiometric_calibration).T
   
                # Reverse channels, catch NaNs, and write
                frame[sp.logical_not(sp.isfinite(frame))]=0
                if config.reverse_channels:
                    frame = sp.flip(frame, axis=0)

                # Clip the channels to the appropriate size, if needed
                if args.extract_subframe:
                    frame = frame_extract(frame, columns = True)

                # Write to file
                sp.asarray(frame, dtype=sp.float32).tofile(fout)
                lines = lines + 1
            
                # Read next chunk
                raw = sp.fromfile(fin, count=config.nraw, dtype=sp.int16)


    # Form output metadata strings
    wl = config.wl_full.copy()
    fwhm = config.fwhm_full.copy()

    if config.extract_subframe:
        ncolumns = last_illuminated_column - first_illuminated_column + 1
        nchannels = last_illuminated_row - first_illuminated_row + 1
        clip_rows = np.arange(first_illuminated_row,
                              last_illuminated_row+1,dtype=int)
        wl = wl[clip_rows]
        fwhm = fwhm[clip_rows]
    else:
        nchannels, ncolumns = native_rows, native_columns

    band_names_string = ','.join(['channel_'+str(i) \
       for i in range(len(wl))])
    fwhm_string =  ','.join([str(w) for w in fwhm])
    wavelength_string = ','.join([str(w) for w in wl])
    

    params = {'lines': lines}
    params.update(**locals())
    with open(config.output_header,'w') as fout:
        fout.write(header_template.format(**params))

    logging.info('Done')


if __name__ == '__main__':

    main()
