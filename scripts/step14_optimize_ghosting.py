# David R Thompson
import os
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
if os.path.exists('/Users/'):
    basedir = '/Users/drt/data/21EMIT/20211210_ghosts/optimization/'
    exe = '/Users/drt/src/emit-sds-l1b/utils/optimizeghost.py'
    cmd = 'bash '+templatefile
    datadir = '/Users/drt/src/emit-sds-l1b/data/'
else:
    basedir = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/'
    exe = '/home/drt/src/emit-sds-l1b/utils/optimizeghost.py'
    srun_flags = '-N 1 -n 1 -c 40 --mem=180G' 
    cmd = 'sbatch '+srun_flags+' '+templatefile
    datadir = '/home/drt/src/emit-sds-l1b/data/'

test_frame  = basedir+'test_frame'
test_frame2 = basedir+'test_frame2'
test_frame3 = basedir+'test_frame3'
test_frame4 = basedir+'test_frame4'
infile = datadir+'emit_ghost.json'
outfile = datadir+'emit_ghost_optimized.json'
cmd = 'python %s %s %s %s'%(exe,infile,test_frame3,outfile)

template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
print(cmd)
os.system(cmd)

