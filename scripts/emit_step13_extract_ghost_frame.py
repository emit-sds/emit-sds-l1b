# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from glob import glob
import sys, os, os.path

indir = '/beegfs/scratch/drt/20220203_EMIT_TVAC4/ghost/' #/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/'
outdir = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/'
infiles = glob(indir+'*pedestal.hdr')#['20211117_084042_UTC_LinearitySphere_Field240_Step15p2mm_PD1000p0candelam2_clip_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr','20211117_044932_UTC_LinearitySphere_Field840_Step14p4mm_PD1616p0candelam2_clip_darksub_pedestal_badfix_osffix_linear_scatterfix.hdr']


for file_index, infilepath in enumerate(infiles):
  infile = os.path.split(infilepath)[-1]
  outfile = 'test_frame_%i.hdr' % (file_index)
  I = envi.open(indir + infile)
  x = I.load()
  frame = np.squeeze(x.mean(axis=0))
  frame = frame.T
  frame = np.array(frame,dtype=np.float32)
  envi.save_image(outdir+outfile,frame,ext='',force=True)



 
