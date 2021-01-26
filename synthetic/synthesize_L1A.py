#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import json
import logging
import argparse
import random
from numpy.random import normal
from numpy.random import randint
from scipy.interpolate import interp1d

linearity_header_string = """ENVI
description = {{Linearity Correction File}}
samples = 65536
lines = 1
bands = 0
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

bad_header_string = """ENVI
description = {{Bad Pixel Map}}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

dark_header_string = """ENVI
description = {{Dark Frame}}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

flat_header_string = """ENVI
description = {{Flat Field}}
samples = {rows}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

crf_header_string = """ENVI
description = {{CRF Correction}}
samples = {cols}
lines = {cols}
bands = 1
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""

srf_header_string = """ENVI
description = {{SRF Correction}}
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

    with open(args.config_file,'r') as fin:
        config = json.load(fin)

    for q in [config['radiometric_coefficient_file'],
              config['spectral_calibration_file'],
              config['flat_field_file'],
              config['crf_correction_file'],
              config['srf_correction_file'],
              config['bad_element_file'],
              config['linearity_file'],
              config['dark_frame_file'],
              config['output_raw_file']]:
        if os.path.exists(q):
             print(q+' exists')
             #sys.exit(1)

    # Make Radiometric Coefficients
    masked = np.array([int(c) for c in config['dark_channels']])
    rows, cols = 328, 1280
    dispersion = 7.4296517
    wl = 265+np.arange(328)*dispersion
    fwhm = np.ones(rows) * (wl[1] - wl[0]) * 1.1
    rccs = np.ones(rows) * 0.01
    rccs[masked]=0.0
    uncerts = np.ones(rows) * 0.000001
    chn = np.arange(rows)

    np.savetxt(config['radiometric_coefficient_file'],np.c_[rccs, uncerts, chn])
    np.savetxt(config['spectral_calibration_file'],np.c_[chn, wl/1000.0, fwhm/1000.0])

    bad = np.zeros((rows,cols), dtype=np.float32)
    nbad = 6
    for n in range(nbad):
        bad[randint(rows),randint(cols)] = -1
    np.array(bad,dtype=np.uint16).tofile(config['bad_element_file'])
    with open(config['bad_element_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(bad_header_string.format(**lcl))

    dark_noise = 10
    dark = np.array(normal(0, dark_noise, (rows,cols)), dtype=np.float32)
    dark = np.tile(dark.reshape((1,rows,cols)),(2,1,1))
    dark.tofile(config['dark_frame_file'])
    with open(config['dark_frame_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(dark_header_string.format(**lcl))
        
    srf = np.array(np.eye(rows), dtype=np.float32)
    srf.tofile(config['srf_correction_file'])
    with open(config['srf_correction_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(srf_header_string.format(**lcl))
        
    crf = np.array(np.eye(cols), dtype=np.float32)
    crf.tofile(config['crf_correction_file'])
    with open(config['crf_correction_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(crf_header_string.format(**lcl))
 
    flat_noise = 0.01
    flat = np.array(normal(1, flat_noise, (rows,cols)), dtype=np.float32)
    flat = np.tile(flat.reshape((1,rows,cols)),(2,1,1))
    flat[1,:,:] = 0.001 # rough uncertainty
    flat.tofile(config['flat_field_file'])
    with open(config['flat_field_file']+'.hdr','w') as fout:
        lcl = locals()
        fout.write(flat_header_string.format(**lcl))
            
    count = np.arange(2**16, dtype=np.uint16).tofile(config['linearity_file'])
    with open(config['linearity_file']+'.hdr','w') as fout:
        fout.write(linearity_header_string) 

    obs_file = envi.open(config['input_obs_file']+'.hdr')
    obs_metadata = obs_file.metadata
    obs_metadata['samples'] = cols

    loc_file = envi.open(config['input_loc_file']+'.hdr')
    loc_metadata = loc_file.metadata
    loc_metadata['samples'] = cols

    in_obs = int(obs_metadata['bands'])
    in_loc = 3
  
    rdn_file = envi.open(config['input_rdn_file']+'.hdr')
    rdn_metadata = rdn_file.metadata.copy()
    in_wl = np.array([float(w) for w in rdn_metadata['wavelength']])
    in_lines = int(rdn_metadata['lines'])
    in_samples = int(rdn_metadata['samples'])
    in_bands = int(rdn_metadata['bands'])

    raw_metadata = rdn_metadata.copy()
    raw_metadata['data type'] = 2 
    raw_metadata['samples'] = cols
    raw_metadata['bands'] = rows

    if any([not q.endswith('.img') for q in [config['output_raw_file'],
        config['output_obs_file'],config['output_loc_file']]]):
          raise ValueError('Filenames must end in .img')

    envi.write_envi_header(config['output_raw_file'].replace('.img','.hdr'),
        raw_metadata)
    envi.write_envi_header(config['output_obs_file'].replace('.img','.hdr'),
        obs_metadata)
    envi.write_envi_header(config['output_loc_file'].replace('.img','.hdr'),
        loc_metadata)
    
    with open(config['output_raw_file'],'wb') as raw_out:
      with open(config['output_obs_file'],'wb') as obs_out:
        with open(config['output_loc_file'],'wb') as loc_out:
          with open(config['input_rdn_file'],'rb') as rdn_in:
            with open(config['input_obs_file'],'rb') as obs_in:
              with open(config['input_loc_file'],'rb') as loc_in:

                for line in range(in_lines):

                    if line%100==0:
                        print(line+1,'/',in_lines)
                    rdn = np.fromfile(rdn_in, count=in_samples*in_bands, dtype=np.float32)
                    rdn = rdn.reshape((in_bands, in_samples)) 
                    rdn_resamp = np.zeros((rows, in_samples))
                    for i in range(in_samples):
                       rdn_resamp[:,i] = interp1d(in_wl, rdn[:,i], bounds_error=False,
                           fill_value=0)(wl)

                    # Thanks to panda-34 of stackoverflow.com
                    desired_shape = np.array((rows, cols))
                    extant_shape = np.array((rdn_resamp.shape[0], rdn_resamp.shape[1]))
                    reps = int(desired_shape[1]/extant_shape[1])+1
                    rdn_resamp = np.tile(rdn_resamp,[1, reps])[:,:desired_shape[1]]
                    raw = np.array(rdn_resamp/(flat[0,:,:].T*rccs).T + dark[0,:,:], 
                        dtype=np.int16)
                    if config['reverse_channels']:
                       raw = np.flipud(raw)
                    raw.tofile(raw_out)
                    
                    obs = np.fromfile(obs_in, count=in_samples*in_obs, dtype=np.float64)
                    obs = obs.reshape((in_samples, in_obs))
                    obs = np.tile(obs,[reps,1])[:desired_shape[1],:]
                    obs.tofile(obs_out)
                    
                    loc = np.fromfile(loc_in, count=in_samples*3).reshape((in_samples, 3))
                    loc = loc.reshape((in_samples, in_loc))
                    loc = np.tile(loc,[reps,1])[:desired_shape[1],:]
                    loc[:,0] = np.linspace(loc[0,0], loc[-1,0], cols)
                    loc[:,1] = np.linspace(loc[0,1], loc[-1,1], cols)
                    loc.tofile(loc_out)

    print('done') 

if __name__ == '__main__':

    main()
