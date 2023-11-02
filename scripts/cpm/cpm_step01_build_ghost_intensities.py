# David R Thompson
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from astropy import modeling
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter

basedir = '/beegfs/scratch/drt/20230901_CPM_GhostStray/20230901_113900_UTC_GhostMapping_DoubleMonochromator/'
margin=10


def find_peak(x):
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)
    fitted_model = fitter(model, np.arange(len(x)), x)
    return fitted_model.mean[0]#, fitted_model.amplitude[0], fitted_model.stddev[0]


with open('../../data/cpm/ghost_pointwise.txt','w') as fout:
  for infile in sorted(glob(basedir+'*pedestal.hdr')):
    print(infile)
    I = envi.open(infile).load()
    ghost = np.squeeze(I[:,465,:])
    source = np.squeeze(I[:,180,:])
    #sourcechan = np.argmax(np.mean(source,axis=0))
    sourcechan = find_peak(np.mean(source,axis=0))
    sourceidx = int(round(sourcechan))
    use = source[:,sourceidx] > 1000
    if sum(use)<5:
        continue
    source_series = np.mean(source[use,:],axis=0)[(sourceidx-margin):(sourceidx+margin+1)]
    done = False
    while True:
        ghostchan = find_peak(np.mean(ghost[use,:],axis=0))
        ghostidx = int(round(ghostchan))
        #ghostchan = np.argmax(np.mean(ghost[use,:],axis=0))
        ghost_series = np.mean(ghost[use,:],axis=0)[(ghostidx-margin):(ghostidx+margin+1)]
        if sum(ghost_series)<7:
        #if sum(ghost_series)<10:
            break
        ratio = sum(ghost_series) / sum(source_series)
        fout.write('%10.8f %10.8f %10.8f %10.8f\n'%\
            (sourcechan, ghostchan, sum(ghost_series), ratio))
        ghost[:,(ghostidx-margin):(ghostidx+margin+1)] = 0
       #if max(ghost_series)>5:
       #  plt.plot(ghost_series)
       #  plt.show()

