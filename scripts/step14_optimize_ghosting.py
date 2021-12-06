# David R Thompson
import os
import numpy as np
import pylab as plt
from spectral.io import envi

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
infile = '/home/drt/src/emit-sds-l1b/data/emit_ghost.json'
outfile = '/home/drt/src/emit-sds-l1b/data/emit_ghost_optimized.json'
#cmd = 'python /home/drt/src/emit-sds-l1b/utils/optimizeghost.py %s %s %s'%(infile,test_frame,outfile)
outfile2 = '/home/drt/src/emit-sds-l1b/data/emit_ghost_optimized_expressive.json'
cmd = 'python /home/drt/src/emit-sds-l1b/utils/optimizeghost.py --expressive %s %s %s %s %s'%(outfile,test_frame,test_frame2,test_frame3,outfile2)
template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
srun_flags = '-N 1 -n 1 -c 40 --mem=180G'# --partition=patient' 
cmd = 'sbatch '+srun_flags+' '+templatefile
print(cmd)
os.system(cmd)

