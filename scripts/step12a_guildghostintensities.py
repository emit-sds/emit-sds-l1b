# David R Thompson
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi

basedir = '/beegfs/scratch/drt/20211114_EMIT_Infield/20211115_InFieldScatter/'
margin=10

for infile in sorted(glob(basedir+'*badfix.hdr')):
  I = envi.open(infile).load()
  ghost = np.squeeze(I[:,358,:])
  source = np.squeeze(I[:,940,:])
  sourcechan = np.argmax(np.mean(source,axis=0))
  use = source[:,sourcechan] > 1000
  if sum(use)<10:
      continue
  source_series = np.mean(source[use,:],axis=0)[(sourcechan-margin):(sourcechan+margin+1)]
  done = False
  while True:
      ghostchan = np.argmax(np.mean(ghost[use,:],axis=0))
      ghost_series = np.mean(ghost[use,:],axis=0)[(ghostchan-margin):(ghostchan+margin+1)]
      if sum(ghost_series)<10:
          break
      ratio = sum(ghost_series) / sum(source_series)
      print(sourcechan, ghostchan, sum(ghost_series), ratio)
      ghost[:,(ghostchan-margin):(ghostchan+margin+1)] = 0
     #if max(ghost_series)>5:
     #  plt.plot(ghost_series)
     #  plt.show()

