# David R Thompson
import sys, os
import numpy as np

# collect linearity curves
cmd = 'python ~/src/emit-sds-l1b/utils/combinelinearity.py '

if False:
  for fieldpoint in [15,90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]:
    fp = str(fieldpoint)
    cmd = cmd + '/beegfs/scratch/drt/20211115_EMIT_Linearity/Field_'+fp+'_curves '
  cmd = cmd + '../data/EMIT_LinearityBasis_20211118'
  print(cmd)
  os.system(cmd)

if True:
  for fieldpoint in [165,240,315,390,840,915,990,1065,1140]:
    fp = str(fieldpoint)
    cmd = cmd + '/beegfs/scratch/drt/20211115_EMIT_Linearity/Field_'+fp+'_fullrangecurves '
  cmd = cmd + '../data/EMIT_LinearityBasis_20211215'
  print(cmd)
  os.system(cmd)
