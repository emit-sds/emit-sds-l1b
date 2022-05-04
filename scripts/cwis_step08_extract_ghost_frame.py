# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi


indir = '/beegfs/scratch/drt/20220112_CWIS2/20220127_IntSph_Ghosts/'
outdir = '/beegfs/scratch/drt/20220112_CWIS2/optimization/'

infiles = ['20220127_IntSp10500_FilterRed14_1mmKG4',    
          '20220127_IntSp10500_FilterRed260_1mmKG4',  
          '20220127_IntSp10700_FilterBrown376',       
          '20220127_IntSp10700_FilterRed36_1mmKG4',   
          '20220127_IntSp12000_Filter1mmBG14_1mmKG4', 
          '20220127_IntSp12000_Filter2mmBG26_1mmKG4', 
          '20220127_IntSp12200_FilterLP514_1mmKG4',   
          '20220127_IntSp12600_FilterLP1319_1mmKG4',  
          '20220127_IntSp13000_Filter1mmBG34_1mmKG4',  
          '20220127_IntSp16000_Filter2mmBG34_1mmKG4',  
          '20220127_IntSp18000_FilterBrown376_1mmKG4', 
          '20220127_IntSp19500_FilterGaAs_1mmKG4',     
          '20220127_IntSp19500_FilterSilicon_1mmKG4',  
          '20220127_IntSp6000_FilterRed260',           
          '20220127_IntSp6700_FilterRed14',            
          '20220127_IntSph19900_Filter2mmKG2',         
          '20220127_IntSph19900_Filter3mmKG4',         
          '20220127_IntSph19900_Filter59044',
          '20220127_IntSph5300',
          '20220127_IntSph6200_FilterRed36',
          '20220127_IntSph6500_Filter1mmBG14',
          '20220127_IntSph6500_Filter2mmBG26',
          '20220127_IntSph7300_Filter1mmBG34',
          '20220127_IntSph8700_Filter2mmBG34',
          '20220127_IntSph9500_Filter1mmKG4']

infiles = []
for file_index, infile in enumerate(infiles):
  outfile = 'test_frame_%s.hdr' % (infile)
  I = envi.open(indir + infile + '_darksub_pedestal.hdr')
  x = I.load()
  lines, samples, bands = x.shape
  frame = np.squeeze(x.mean(axis=0))
  frame = frame.reshape((1,samples,bands))
  envi.save_image(outdir+outfile,frame,ext='',interleave='bil',force=True)

indir = '/beegfs/scratch/drt/20220112_CWIS2/20211213_crf/'
infile = '20211213_CRF_Test'

outfile = 'test_frame_%s.hdr' % (infile)
I = envi.open(indir + infile + '_darksub_pedestal.hdr')
x = I.load()
lines, samples, bands = x.shape
frame = np.squeeze(x.mean(axis=0))
frame = frame.reshape((1,samples,bands))
frame[:,78:1217] = 0
envi.save_image(outdir+outfile,frame,ext='',interleave='bil',force=True)

 
