# David R Thompson
from glob import glob
import sys, os

for filepath in glob('/beegfs/scratch/drt/20220112_CWIS2/linearity/PD*pedestal'):
  toks = filepath.split('_')
  toks[-3] = 'Field365_PD' + toks[-3] + 'p0candelam2'
  newpath = '_'.join(toks)
  cmd = 'cp %s %s' % (filepath, newpath)
  print(cmd)
  os.system(cmd)
  cmd = 'cp %s.hdr %s.hdr' % (filepath, newpath)
  print(cmd)
  os.system(cmd)
