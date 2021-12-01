# David R Thompson
import sys, os
import numpy as np
from glob import glob
from spectral.io import envi

outframes = []

files = sorted(glob('/beegfs/scratch/drt/20211130_EMIT_Ghost/20211117_194709_UTC_InFieldScatter/*linear'))
for infile in files:
    X = envi.open(infile+'.hdr').load()
    X = np.mean(X,axis=0)
    print(X.shape)
    outframes.append(X)

files = sorted(glob('/beegfs/scratch/drt/20211130_EMIT_Ghost/20211116_045752_UTC_InFieldScatter/*linear'))
for infile in files:
    freq = float(infile.split('/')[-1].split('_')[4].replace('nm','').replace('p','.'))
    print(infile,freq)
    if freq<=490:
        continue
    X = envi.open(infile+'.hdr').load()
    X = np.mean(X,axis=0)
    outframes.append(X)

X = np.array(outframes, dtype=np.float32)
envi.save_image('/beegfs/scratch/drt/20211130_EMIT_Ghost/combined_frames.hdr',X,interleave='bil',ext='',force=True)

