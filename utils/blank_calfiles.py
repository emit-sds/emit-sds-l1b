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
import random
from numpy.random import normal


linearity_header_string = """ENVI
description = {Linearity Correction File}
samples = 65536
lines = 1
bands = 0
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

bad_header_string = """ENVI
description = {Bad Pixel Map}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

dark_header_string = """ENVI
description = {Dark Frame}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

flat_header_string = """ENVI
description = {Flat Field}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

crf_header_string = """ENVI
description = {CRF Correction}
samples = {cols}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

srf_header_string = """ENVI
description = {SRF Correction}
samples = {rows}
lines = {rows}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""




def main():

    description = "Synthesize calibration files and DNs"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config_file')
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
    rows, cols = 328, 1280
    dispersion = 7.4296517
    wl = 265+np.arange(328)*dispersion
    fwhm = (wl[1] - wl[0]) * 1.1
    rccs = np.array(rows) * 0.0001
    uncerts = np.ones(rows) * 0.000001
    chn = np.arange(rows)

    np.savetxt(np.c_[rccs, uncerts, chn], 
        config['radiometric_coefficient_file'])

    np.savetxt(np.c_[chn, wl,fwhm], 
        config['spectral_calibration_file'])

    bad = np.zeros((rows,cols), dtype=np.float32)
    nbad = 6
    for n in nbad:
        bad[random.randint(rows),random.randint(cols))] = -1
    np.tofile(config['bad_element_file'], bad)
    with open(config['bad_element_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(bad_header_string.format(**lcl))

    dark_noise = 10
    dark = np.array(normal(0, dark_noise, (rows,cols)), dtype=np.float32)
    np.tofile(config['dark_frame_file'], dark)
    with open(config['dark_frame_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(dark_header_string.format(**lcl))
        
    srf = np.array(np.eye(rows), dtype=np.float32)
    np.tofile(config['srf_correction_file'], dark)
    with open(config['srf_correction_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(srf_header_string.format(**lcl))
        
    crf = np.array(np.eye(rows), dtype=np.float32)
    np.tofile(config['crf_correction_file'], dark)
    with open(config['crf_correction_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(crf_header_string.format(**lcl))
 
    flat_noise = 0.01
    flat = np.array(normal(1, flat_noise, (rows,cols)), dtype=np.float32)
    np.tofile(config['flat_field_file'], flat)
    with open(config['flat_field_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(flat_header_string.format(**lcl))
            
    count = sp.arange(2**16, dtype=sp.uint16).tofile(args['linearity_file'])
    with open(config['linearity_file'],'w') as fout:
        fout.write(header_string) 

    obs_file = envi.open(args.obs_source_file+'.hdr')
    obs_metadata = obs_file.meta()
    in_obs = obs_metadata['bands']
  
    rdn_file = envi.open(args.rdn_source_file+'.hdr')
    rdn_metadata = rdn_file.meta()
    raw_metadata = rdn_metadata.copy()
    in_lines = rdn_metadata['lines']
    in_samples = rdn_metadata['samples']
    in_bands = rdn_metadata['bands']
    raw_metadata['data type'] = 2 
    raw_metadata['samples'] = cols
    envi.write_header(config['input_file']+'.hdr',raw_metadata,ext='',force=False)
    
    with open(config['simulated_raw_file'],'wb') as rdn_out:
      with open(config['simulated_obs_file'],'wb') as obs_out:
        with open(config['simulated_loc_file'],'wb') as loc_out:
          with open(config['input_rdn_file'],'rb') as rdn_in:
            with open(config['input_obs_file'],'rb') as obs_in:
              with open(config['input_loc_file'],'rb') as loc_in:

                rdn = np.fromfile(rdn_in, count=in_samples*in_bands, dtype=np.float32)
                rdn = rdn.reshape((in_bands, in_samples)) 
                rdn_resamp = np.zeros((rows, in_samples))
                for i in range(len(rdn_resamp)):
                   rdn_resamp[:,i] = interp1d(in_wl, rdn, bounds_error=False,
                       fill_value='extrapolate')
                valid = np.where(np.all(rdn>-9990,axis1))[0]
                # Thanks to panda-34 of stackoverflow.com
                rdn_valid = rdn_resamp[:,valid:]
                desired_shape = np.array((rows, cols))
                pads = tuple((0, i) for i in (desired_shape-rdn_valid.shape))
                rdn_out = np.pad(rdn_valid, pads, mode="wrap")

                obs = np.fromfile(obs_in, count=in_samples*in_obs, dtype=np.float32)
                obs = obs.reshape((in_samples, in_obs))
                obs_valid = obs[valid:,:]
                desired_obs_shape = np.array((cols,in_obs))
                pads = tuple((0, i) for i in (desired_obs_shape-obs_valid.shape))
                obs_out = np.pad(obs_valid, pads, mode="wrap")

                loc = np.fromfile(loc_in, count=in_samples*3).reshape((in_samps,out_samps))
                loc = loc.reshape((in_samples, in_loc))
                loc_valid = loc[valid:,:]
                desired_loc_shape = np.array((cols, 3))
                pads = tuple((0, i) for i in (desired_loc_shape-loc_valid.shape))
                loc_out = np.pad(loc_valid, pads, mode="wrap")
                loc_out[:,0] = np.linspace(loc_out[0,0], loc_out[-1,0], cols)
                loc_out[:,1] = np.linspace(loc_out[0,1], loc_out[-1,1], cols)


                obs_out = np.pad(arr, pads
    
    print('done') 

if __name__ == '__main__':

    main()
