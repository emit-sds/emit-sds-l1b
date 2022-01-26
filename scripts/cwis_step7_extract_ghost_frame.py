# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi


indir = '/beegfs/scratch/drt/20220112_CWIS2/ghost_measurements/'
outdir = '/beegfs/scratch/drt/20220112_CWIS2/optimization/'
infiles = ['20220119_Filter_BP1570_addSilicon',
           '20220119_Filter_BP1570',
           '20220119_Filter_GaAs',
           '20220119_Filter_LP561_addKG2_9mm',
           '20220119_Filter_MysteryCyan',
           '20220119_Filter_SP785_addKG2_12mm',
           '20220119_Filter_SP785',
           '20220119_Filter_addKG2_4andhalfmm',
           '20220119_Filter_addKG4_3mm',
           '20220119_PD4000']

for file_index, infile in enumerate(infiles):
  outfile = 'test_frame_%i.hdr' % (file_index)
  I = envi.open(indir + infile + '_darksub_pedestal.hdr')
  x = I.load()
  lines, samples, bands = x.shape
  frame = np.squeeze(x.mean(axis=0))
  frame = frame.reshape((1,samples,bands))
  #frame = np.array(frame,dtype=np.float32)
  envi.save_image(outdir+outfile,frame,ext='',interleave='bil',force=True)


