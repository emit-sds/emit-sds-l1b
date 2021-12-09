# David R Thompson
import os
import numpy as np
import pylab as plt
from spectral.io import envi


# segment the ghost orders into parts, based on a manual assessment of discontinuities
# and inflections in the ghost spectrum.  Write the results as a new ghost file
cmd = 'python ../utils/segmentghost.py ../data/emit_ghost.json ../data/emit_ghost_segmented.json'
os.system(cmd)

batch_template='''#!/bin/sh
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1 
ray stop
ps -ef | grep ray | awk '{print $1}' | xargs kill
%s 
'''

templatefile = 'batch_run_delete.sh'
test_frame = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame.hdr'
test_frame2 = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame2.hdr'
test_frame3 = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame3.hdr'
test_frame4 = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/test_frame4.hdr'
infile = '/home/drt/src/emit-sds-l1b/data/emit_ghost_segmented.json'
outfile = '/home/drt/src/emit-sds-l1b/data/emit_ghost_optimized.json'
cmd = 'python /home/drt/src/emit-sds-l1b/utils/optimizeghost.py %s %s %s %s'%(infile,test_frame2,test_frame4,outfile)
template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
srun_flags = '-N 1 -n 1 -c 40 --mem=180G'# --partition=patient' 
cmd = 'sbatch '+srun_flags+' '+templatefile
print(cmd)
os.system(cmd)

