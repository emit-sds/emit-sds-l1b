# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from glob import glob
import sys, os, os.path

indir = '/beegfs/scratch/drt/20230901_CPM_GhostStray/20230901_072600_UTC_GhostMapping_Whitelight/'
outdir = indir+'/optimization/'
infiles = glob(indir+'*pedestal.hdr')


for file_index, infilepath in enumerate(infiles):
  infile = os.path.split(infilepath)[-1]
  outfile = 'test_frame_%i.hdr' % (file_index)
  I = envi.open(indir + infile)
  samples = int(I.metadata['samples'])
  bands = int(I.metadata['bands'])
  x = I.load()
  frame = np.squeeze(x.mean(axis=0))
  frame = frame.reshape((1,samples,bands))
  frame = np.array(frame,dtype=np.float32)
  envi.save_image(outdir+outfile,frame,ext='',interleave='bil',force=True)



 
