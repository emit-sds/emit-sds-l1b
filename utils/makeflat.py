# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from skimage.filters import threshold_otsu
from numpy import nanmedian
import json
from numba import jit
from lowess import lowess
from PIL import Image


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def moving_average(x, w=5):
      return np.convolve(x, np.ones(w), 'same') / w


def polymax(y, plot=False):
    series = moving_average(y)
    ctr = np.argmax(series)
    halfwid = 16
    segment = y[max(0,ctr-halfwid):min(ctr+halfwid+1,len(y)-1)]
    x = np.arange(len(segment))
    p = np.polyfit(x,segment,6)
    noisefree = np.polyval(p,x)
    if plot:
        plt.plot(x, segment, 'ko')
        plt.plot(x, noisefree, 'r')
        plt.show()
    return noisefree.max(), np.std(noisefree-segment)


# Reference columns of the focal plane array used for
# radiometric calibration.  Avoid the center (due to 
# symmetric ghosting) and avoid the divot from 1015-1035.
reference_cols = np.concatenate((np.arange(140,340),
                            np.arange(940,1015),
                            np.arange(1035,1140)),axis=0)



def main():

    description = "Calculate Flat field"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--cue_channel',default=50,type=int)
    parser.add_argument('--background',type=str)
    parser.add_argument('--config',type=str)
    parser.add_argument('--mask_image',type=str,default=None)
    parser.add_argument('output')
    args = parser.parse_args()

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
    meta = {'interleave':'bsq', 'data type':4}
        
    flat = np.ones((rows,columns)) * -9999
    noise = np.ones((rows,columns)) * -9999
    DN_average, DN_noise = [],[]

    if args.mask_image is not None:
 
        mask = np.asarray(Image.open(args.mask_image))
        print(mask.shape) 
        if len(mask.shape)>2:
            mask = mask.max(axis=2)
        if mask.shape[0] != lines or mask.shape[1] != columns:
            raise IndexError('mask does not match image')

        n = max(mask.sum(axis=0))
        foreground = np.zeros((rows,columns))
        foreground_counts = np.zeros((rows,columns))
        foreground_sq = np.zeros((rows,columns))
        background = np.zeros((rows,columns))
        background_counts = np.zeros((rows,columns))
        background_sq = np.zeros((rows,columns))

        with open(args.input,'rb') as fin:
        
            # Accumulate n brightest observations of the source
            for line in range(lines):
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                use = np.where(mask[line,:]>0)
                if len(use)<1:
                    continue
                foreground[:,use] = foreground[:,use] + frame[:,use]
                foreground_sq[:,use] = foreground_sq[:,use] + pow(frame[:,use],2)
                foreground_counts[:,use] = foreground_counts[:,use] + 1
            foreground = foreground / foreground_counts   
            foreground_sq = foreground_sq / foreground_counts   

        if args.background is not None:
            with open(args.background,'rb') as fin:
                for line in range(lines):
                    frame = np.fromfile(fin, count=nframe, dtype=dtype)
                    frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                    use = np.where(mask[line,:]>0)
                    if len(use)<1:
                        continue
                    background[:,use] = background[:,use] + frame[:,use]
                    background_sq[:,use] = background_sq[:,use] + pow(frame[:,use],2)
                    background_counts[:,use] = background_counts[:,use] + 1
                background = background / background_counts   
                background_sq = background_sq / background_counts   

        foreground_sd = np.sqrt(foreground_sq + pow(foreground,2))
        background_sd = np.sqrt(background_sq + pow(background,2))

        for row in range(rows):
        
           flat[row,:] = foreground[row,:] - background[row,:]
           noise[row,:] = np.sqrt(pow(foreground_sd[row,:],2)-pow(background_sd[row,:],2))
        
           ref = nanmedian(flat[row, reference_cols])
           ref_noise = nanmedian(noise[row, reference_cols])
           print('row',row,'reference average is',ref)
           flat[row,:] = ref / flat[row,:] 
           DN_average.append(ref)
           DN_noise.append(ref_noise)

    else:

        foreground = np.ones((lines,rows,columns))
        background = np.ones((lines,rows,columns))
        with open(args.input,'rb') as fin:
        
            # Accumulate n brightest observations of the source
            for line in range(lines):
                frame = np.fromfile(fin, count=nframe, dtype=dtype)
                frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
                foreground[line,:,:] = frame
               
        if args.background is not None:
            with open(args.background,'rb') as fin:
        
                # Accumulate n brightest observations of the background
                for line in range(lines):
                    bg = np.fromfile(fin, count=nframe, dtype=dtype)
                    bg = np.array(bg.reshape((rows, columns)), dtype=np.float32)
                    background[line,:,:] = bg
        
        for row in range(rows):
           for col in range(columns):
        
               y = np.squeeze(foreground[:,row,col])
               bg_y = np.squeeze(background[:,row,col])
               fg, resid_fg = polymax(y,plot=False)#(row==150 and col==200))
               bg, resid_bg = polymax(bg_y,plot=False)#(row==150 and col==200))
               flat[row,col] = fg - bg
               noise[row,col] = resid_fg
        
           ref = nanmedian(flat[row, reference_cols])
           ref_noise = nanmedian(noise[row, reference_cols])
           print('row',row,'reference average is',ref)
           flat[row,:] = ref / flat[row,:] 
           DN_average.append(ref)
           DN_noise.append(ref_noise)

    flat[np.logical_not(np.isfinite(flat))] = -9999
    meta['average_DNs'] = np.array(DN_average)
    meta['stdev_DNs'] = np.array(DN_noise)
    envi.save_image(args.output+'.hdr',np.array(flat,dtype=np.float32),
        metadata=meta,ext='',force=True)



if __name__ == '__main__':

    main()
