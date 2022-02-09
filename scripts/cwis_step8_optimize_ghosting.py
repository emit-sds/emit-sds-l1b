# David R Thompson
import os, os.path
import numpy as np
import pylab as plt
from spectral.io import envi


# segment the ghost orders into parts, based on a manual assessment of discontinuities
# and inflections in the ghost spectrum.  Write the results as a new ghost file
#cmd = 'python ../utils/segmentghost.py ../data/emit_ghost.json ../data/emit_ghost_segmented.json'
#os.system(cmd)

batch_template='''#!/bin/sh
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1 
ray stop
ps -ef | grep ray | awk '{print $1}' | xargs kill
%s 
'''

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


templatefile = 'batch_run_delete.sh'
mydir = os.path.split(os.path.abspath(__file__))[0]+'/'
exe = mydir + '/../utils/optimizeghost.py'

basedir = '/beegfs/scratch/drt/20220112_CWIS2/optimization/'
srun_flags = '-N 1 -n 1 -c 40 --mem=180G' 
cmd = 'sbatch '+srun_flags+' '+templatefile
datadir = mydir+'../data/'

infile = datadir+'cwis_ghost.json'
outfile = datadir+'../data/CWIS_GhostMap_20220125.json'
cmd = 'srun '+srun_flags+' python %s %s '%(exe,infile)

for infile in infiles:
   cmd = cmd +basedir+('test_frame_%s'%infile)+' '
cmd = cmd + outfile
template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
print(cmd)
os.system(cmd)


