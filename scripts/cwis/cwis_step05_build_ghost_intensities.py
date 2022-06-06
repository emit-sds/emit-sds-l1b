# David R Thompson
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from astropy import modeling
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter

basedir = '/beegfs/scratch/drt/20220112_CWIS2/hiss/'
margin=20

infiles = ['20220124t1530avg10_scan550to700nm_SAT_clip_darksub.hdr',
           '20220124t1535avg10_scan700to850nm_SAT_clip_darksub.hdr',
           '20220124t1535avg10_scan850to1200nm_SAT_clip_darksub.hdr',
           '20220124t1455avg10_scan1200to1700nm_clip_darksub.hdr',
           '20220124t1505avg10_scan1700to2000nm_clip_darksub.hdr',
           '20220124t1510avg10_scan2000to2500nm_clip_darksub.hdr']

infiles = [basedir + f for f in infiles]

def find_peak(x):
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)
    fitted_model = fitter(model, np.arange(len(x)), x)
    return fitted_model.mean[0]#, fitted_model.amplitude[0], fitted_model.stddev[0]


with open('../../data/cwis/cwis_ghost_pointwise.txt','w') as fout:
  for infile in infiles:
    print(infile)
    I = envi.open(infile).load()

    orders = None

    for i in range(I.shape[0]):

        ghost = np.squeeze(I[i,931,:])
        source = np.squeeze(I[i,362,:])
        sourcechan = find_peak(source)
        sourceidx = int(round(sourcechan))
        source_series = source[(sourceidx-margin):(sourceidx+margin+1)]

        if sum(source_series)<5:
            continue

       #if sourceidx>199 and sourceidx<213:
       #     continue

        done = False

        while True:

            ghostchan = find_peak(ghost)
            if not np.isfinite(ghostchan):
                break

            ghostidx = int(round(ghostchan))
            ghost_series = ghost[(ghostidx-margin):(ghostidx+margin+1)]

            if len(ghost_series)<20 or max(ghost_series)<20:
                break

           #plt.plot(ghost_series)
           #plt.show()

            ratio = sum(ghost_series) / sum(source_series)
            data = (sourcechan, ghostchan, sum(ghost_series), ratio)

            if orders is None:
                new_order = [data]
                orders = [[data]]
            else:
                for test_order in range(len(orders)-1,-1,-1):
                   lastchan = orders[test_order][-1][1]
                   if abs(ghostchan-lastchan)<10:
                   #if ghostchan > lastchan:
                       orders[test_order].append(data) 
                       break
                   elif test_order == 0:
                       orders.append([data])
            ghost[(ghostidx-margin):(ghostidx+margin+1)] = 0

    if orders is not None:
        for oi, order in enumerate(orders):
            for entry in order:
                outstr = '%10.8f %10.8f %10.8f %10.8f'%entry
                print(oi,outstr)
                fout.write(outstr+'\n')
            fout.write('\n')

