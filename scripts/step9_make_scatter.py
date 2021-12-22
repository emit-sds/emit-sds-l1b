# David R Thompson
import sys, os
import numpy as np
from glob import glob

if False:

    files = glob('/beegfs/scratch/drt/20211114_EMIT_Infield/20211114_InFieldScatter/*linear')
    files.sort()
    cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--target_row','40','--spatial']+files+['>>','spatial_params.txt']
    os.system(' '.join(cmds))


    files = glob('/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/*linear')
    files.sort()
    cmds = ['python','/home/drt/src/emit-sds-l1b/utils/makescatter.py','--target_row','940']+files+['>>','spectral_params.txt']
    os.system(' '.join(cmds))

if True:

    cmd = 'python ../utils/combinescatter.py --spatial  ../scripts/spatial_params.txt ../data/EMIT_SpatialScatter_20211122' 
    os.system(cmd)

    cmd = 'python ../utils/combinescatter.py ../scripts/spectral_params.txt ../data/EMIT_SpectralScatter_20211122' 
    os.system(cmd)
 
 
 
 
 
 
 
 
 
