#!/usr/bin/env python
# David R Thompson 
# Spring 2015
# Jet Propulsion Laboratory, California Institute of Technology

import os, sys, argparse, time
import spectral
from scipy.optimize import minimize
import pylab as plt
from scipy.interpolate import splrep, splev
from scipy import ndimage 
from skimage import filters
import numpy as np
import spectral.io.envi as envi
from numba import jit
from fpa import FPA
from scipy.ndimage import gaussian_filter

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def main():

    description = 'Apply a multiplicative flat field (and optionally, an additive offset)';
    parser = argparse.ArgumentParser()

    # Required 
    parser.add_argument('input', nargs='+', help='Flat field')
    parser.add_argument('--config',type=str)
    parser.add_argument('--highpass',action='store_true')
    parser.add_argument('output_all', help='Output flat field list')
    parser.add_argument('output', help='Output final flat field')
    args = parser.parse_args()
    fpa = FPA(args.config)

    reference_cols = []
    for extrema in fpa.reference_cols:
      reference_cols.extend(np.arange(extrema[0],extrema[1]))
    reference_cols = np.array(reference_cols, dtype=int)

    # Define local variables
    flats = []
    for infile in args.input:
        print(infile)
        inhdr  = find_header(infile)
        img = envi.open(inhdr).load()
        if np.any(img<-9990):
           continue
        band = np.squeeze(img[:,:,50])
        nbands = img.shape[2]
        nsamps = img.shape[1]
       #plt.imshow(band)
       #plt.show()
        edges = abs(ndimage.sobel(ndimage.gaussian_filter(band, 3)))
       #plt.imshow(edges)
       #plt.show()
        thresh = filters.threshold_otsu(edges) 
        thresh = np.percentile(edges,70)
        edge = edges>thresh
       #plt.imshow(edge)
       #plt.show()
        edge = ndimage.binary_dilation(edge)
        bright = np.any(img>35,axis=2)
        use = np.logical_or(edge, bright)
        for i in range(img.shape[0]):
            img[i,use[i,:],:] = np.nan

        flat = np.nanmedian(img,axis=0).T # convert to row, column
        new = flat.copy()
        if args.highpass:
           blur = gaussian_filter(new,(0.4,2))
           new = new / blur
        else:
           for j in range(new.shape[0]):
               knots = np.array([40,330,620,910,1200])
               x = np.arange(new.shape[1])
               sp = splrep(x,new[j,:],t=knots)
              #plt.plot(x,new[j,:])
               new[j,:] = new[j,:]/splev(x,sp) 
              #plt.plot(x,new[j,:])
              #plt.show()
        flat = new
        for row in range(flat.shape[0]):
           ref = np.nanmedian(flat[row, reference_cols])
           flat[row,:] = ref / flat[row,:] 
        flat[np.logical_not(np.isfinite(flat))]= 1
        envi.save_image(infile+'_flat.hdr',np.array(flat,dtype=np.float32),ext='',force=True)
        print(flat)
        flats.append(flat)

    flats = np.array(flats)
    avgflat = np.nanmedian(np.array(flats),axis=0)
    envi.save_image(args.output_all+'.hdr',np.array(flats,dtype=np.float32),ext='',force=True)

    I = envi.open(inhdr)
    meta = I.metadata.copy()
    outhdr = args.output + '.hdr'
    meta['lines'] = nbands
    meta['samples'] = nsamps
    meta['bands'] = 1
    meta['interleave'] = 'bsq'
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")
    with open(args.output,'wb') as fout:
        np.array(avgflat, dtype=np.float32).tofile(fout)
                   
if __name__ == '__main__':
    main()


