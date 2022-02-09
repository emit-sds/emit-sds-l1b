# David R Thompson
import os, os.path
from glob import glob
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

basedir = '/beegfs/scratch/drt/20211130_EMIT_Ghost/optimization/'
srun_flags = '-N 1 -n 1 -c 40 --mem=180G' 
cmd = 'sbatch '+srun_flags+' '+templatefile
datadir = mydir+'../data/'

infiles = [s.replace('.hdr','') for s in glob(basedir+'test_frame_*.hdr')]
infile = datadir+'emit_ghost.json'
outfile = datadir+'../data/EMIT_GhostMap_20220117.json'
cmd = 'srun '+srun_flags+' python %s %s '%(exe,infile)

for infile in infiles:
   cmd = cmd + ' '+infile
cmd = cmd + ' ' + outfile
template = batch_template % cmd
with open(templatefile,'w') as fout:
   fout.write(template)
print(cmd)
os.system(cmd)



