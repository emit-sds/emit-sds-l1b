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

templatefile = 'batch_run_delete.sh'
mydir = os.path.split(os.path.abspath(__file__))[0]+'/'
exe = mydir + '/../utils/optimizeghost.py'

if os.path.exists('/Users/'):
    basedir = '/Users/drt/data/21EMIT/20211210_ghosts/optimization/'
    cmd = 'bash '+templatefile
    datadir = mydir+'../data/'
else:
    basedir = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/'
    srun_flags = '-N 1 -n 1 -c 40 --mem=180G' 
    cmd = 'sbatch '+srun_flags+' '+templatefile
    datadir = mydir+'../data/'

test_frame_0 = basedir+'test_frame_0'
test_frame_1 = basedir+'test_frame_1'

infile = datadir+'emit_ghost.json'
outfile = datadir+'../data/EMIT_GhostMap_20220117.json'
cmd = 'python %s %s %s %s %s'%(exe,infile,test_frame_0,test_frame_1,outfile)
template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
print(cmd)
os.system(cmd)


