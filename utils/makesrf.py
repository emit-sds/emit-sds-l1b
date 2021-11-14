# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from scipy.signal import deconvolve
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from scipy.interpolate import interp1d
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def find_peak(x):
    
    fitter = modeling.fitting.LevMarLSQFitter()

    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)   # depending on the data you need to give some initial values
    fitted_model = fitter(model, np.arange(len(x)), x)
    return fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]


def main():

    description = "Calculate SRFs"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--wavelengths',type=str,default=None)
    parser.add_argument('--monochromator_bandwidth_nm',type=float,default=1.0)
    parser.add_argument('--deconvolve',action='store_true')
    parser.add_argument('--target_index',type=int,default=-1)
    args = parser.parse_args()

    if args.wavelengths is None:
       script = os.path.realpath(__file__)
       directory = os.path.split(script)[0]
       q,wl,fwhm = np.loadtxt(directory+'/../data/EMIT_Wavelengths_20211104.txt').T * 1000.0
    else:
       q,wl,fwhm = np.loadtxt(args.wavelengths).T * 1000.0

    infile = envi.open(find_header(args.input))
 
    if int(infile.metadata['data type']) == 2:
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
    margin=2

    allctrs,alllines = [],[]
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            if args.target_index < 0:
                maxind = np.argmax(np.sum(frame[:330,:],axis=0))
            else:
                maxind = args.target_index
            ctr,amplitude,std = find_peak(frame[:,maxind])
            fwhm = std * 2.0 * np.sqrt(2.0*np.log(2))
            allctrs.append(ctr)
            alllines.append(line)
            #print(line,ctr,fwhm)
   
    chans_per_frame, offset = np.polyfit(alllines,allctrs,1)
    robust_model = RANSACRegressor()
    robust_model.fit(np.array(alllines)[:,np.newaxis],np.array(allctrs))
    chans_per_frame_robust = robust_model.estimator_.coef_[0]
    chans_per_frame = abs(chans_per_frame)
    chans_per_frame_robust = abs(chans_per_frame_robust)
   #print('monochromator velocity:',chans_per_frame,'channels per frame')
   #print('robust monochromator velocity:',chans_per_frame_robust,'channels per frame')

    chans = np.unique([int(round(c)) for c in allctrs])
    sequences = {c:[] for c in chans}
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            maxind = np.argmax(np.sum(frame[:330,:],axis=0))
            if args.target_index < 0:
                maxind = np.argmax(np.sum(frame[:330,:],axis=0))
            else:
                maxind = args.target_index
            for ctr in chans:
                sequences[ctr].append(frame[ctr,maxind])
        
    for ctr, sequence in sequences.items():
         nm_per_channel = abs(wl[ctr] - wl[ctr+1])
         if args.deconvolve:
             # resample to 0.01 nm and deconvolve monochromator
             spacing_frames = 0.01 / nm_per_channel / chans_per_frame_robust
             new_grid = np.linspace(0,len(sequence),  len(sequence)/spacing_frames+1)
             resampled_0p01nm = interp1d(np.arange(len(sequence)), sequence,
                 fill_value='extrapolate',bounds_error=False)(new_grid)
             monochromator = np.ones(int(args.monochromator_bandwidth_nm/0.01))
             monochromator = monochromator / monochromator.sum()
             deconvolved = deconvolve(resampled_0p01nm, monochromator)
             c,amp,std = find_peak(deconvolved[0])
             fwhm = std * 2.0 * np.sqrt(2.0*np.log(2)) # FWHM in grid points
             fwhm = fwhm * 0.01  # FWHM in nm
             fwhm = fwhm / nm_per_channel  # FWHM in channels
         else:
             c,amp,std = find_peak(sequence)
             fwhm = std * 2.0 * np.sqrt(2.0*np.log(2))
             fwhm = fwhm * chans_per_frame_robust # FWHM in channels
             fwhm = fwhm * nm_per_channel # FWHM in nm
         print(ctr,fwhm)
    print('done') 

if __name__ == '__main__':

    main()
