# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from scipy.interpolate import interp1d
from sklearn.linear_model import RANSACRegressor
from sklearn.decomposition import PCA
from numpy import nanmedian
import json
from fpa import FPA
from lowess import lowess


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')




def main():

    description = "Calculate Bad Pixels from average frame"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--threshold',type=float,default=40)
    parser.add_argument('--config',type=str,default=None)
    parser.add_argument('output')
    args = parser.parse_args()
    fpa = FPA(args.config)
    use = np.arange(fpa.first_illuminated_column,fpa.last_illuminated_column+1)

    infile = envi.open(find_header(args.input))
    frame = infile.load()
    std = np.squeeze(frame[:,:,1])
    frame = np.squeeze(frame[:,:,0])
    rows, columns = frame.shape
    mask = np.zeros((rows,columns))
    
    frame = (frame.T - frame[:,use].mean(axis=1)).T
    for col in use:
        spectrum = frame[:,col]
        spectrum[175:184] = interp1d([175,183],[spectrum[175],spectrum[183]])(np.arange(175,184))
        chans = np.arange(rows)
        sm = lowess(spectrum, chans, frac=0.2, return_sorted=False)
        spectrum = spectrum - sm
        #plt.plot(chans,spectrum)
        bad = abs(spectrum)>args.threshold
        #plt.plot(chans[bad],spectrum[bad],'ko')
        #plt.show()      
        mask[bad,col] = 1
        print(sum(bad),' bad pixels in column ',col)

    bads = 0
    bad_map = mask.copy()
    bad_map = np.array(bad_map,dtype=np.int16)
    for column in range(bad_map.shape[1]):
        state_machine = 0
        for row in range(bad_map.shape[0]):
            if mask[row,column]:
                state_machine = state_machine + 1
                bad_map[row,column] = -state_machine
                print(row,column,state_machine)
                bads = bads + 1
            else:
                state_machine = 0
    
    print('total bads:',bads)
    bad_map = bad_map.reshape((rows,columns,1))

    envi.save_image(args.output+'.hdr',
        bad_map, interleave='bsq', ext='', force=True) 


if __name__ == '__main__':

    main()
