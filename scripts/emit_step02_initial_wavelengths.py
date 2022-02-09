# David R Thompson
import pylab as plt
import numpy as np
from scipy.stats import norm
from spectral.io import envi
from scipy.interpolate import interp1d
import astropy.modeling as modeling
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.optimize import minimize
from numpy.random import randn
import json


def find_peak(x, plot=False):
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)   # depending on the data you need to give some initial values
    fitted_model = fitter(model, np.arange(len(x)), x)
    if plot:
         print(fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0])
         plt.plot(x)
         plt.plot(fitted_model(np.arange(len(x))))
         plt.show()
    return fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]


basedir = '/beegfs/scratch/drt/20220124_EMIT_LaserSphere/'
I = envi.open(basedir+'20211112_211905_UTC_LaS_Fields-40-1319_clip_darksub_pedestal_badfix_osffix_linear_scatterfix_ghostfix.hdr')
I = I.load()

# Laser info from 20211112_211905_UTC_LaS_Fields-40-1319
# Wavelengths from Christine.  We start with an initial first guess of the channels
# based on her mesurements.
wavelengths = np.array([1949.10,1064.38,632.8,532.19,405.56]) 
sampling = 8.6
uncert = sampling *  np.array([0.073, 0.007, 0.026, 0.079, 0.077])
channels = np.array([100.36,219.29,277.15,290.62,307.66])
nlasers = len(channels)

# Change to refractive wavelength of vacuum
index_of_refraction = np.array([1.000268,1.000269,1.000271,1.000273,1.000277])
wavelengths = wavelengths * index_of_refraction
 
# These channels are reported in whole-FPA format, but the EMIT FPA begins reading at row 6
row_offset = 6
nrows, ncols = 328, 1280
channels = channels - row_offset 

margin = 5
observed = [[] for c in channels]

for line in range(I.shape[0]):
   frame = np.squeeze(I[line,:,:]).T 
   col,amp,_ = find_peak(frame.mean(axis=0))
   print(col,amp)
   if amp<100:
       continue
   for i, w, chn in zip(range(nlasers), wavelengths, channels):
       idx = np.arange(int(chn-margin),int(chn+margin+1), dtype=int)
       row,_,_ = find_peak(frame[idx,int(round(col))])
       row = row+idx[0] 
       observed[i].append([col,row,w])

x = np.zeros((nlasers,ncols))
y = np.zeros((nlasers,ncols))
for i in range(len(observed)):
    D = np.array(observed[i])
    p = np.polyfit(D[:,0],D[:,1],2)
    x[i,:] = np.polyval(p,np.arange(ncols))
    y[i,:] = wavelengths[i]
    plt.plot(D[:,0],D[:,1]-p[-1],'.')
    plt.plot(np.arange(ncols),x[i,:]-p[-1],'k')
plt.show()


# Perform the fit for each column
ctrs = np.zeros((nrows, ncols, 2))
for i in range(ncols):

  p = np.polyfit(x[:,i],y[:,i],1)
  ctrs[:,i,0] = np.polyval(p,np.arange(nrows))

# Now simulate wavelength error due to uncertain wavelengths
errs = []
ref = 640
for trial in range(1000):
    a = x[:,ref]
    y = wavelengths
    p = np.polyfit(a,y,1)
    y2 = y + randn() * uncert
    p2 = np.polyfit(a,y2,1)
    err = np.polyval(p,np.arange(nrows)) - np.polyval(p2,np.arange(nrows))
    errs.append(err)
errs = abs(np.array(errs)).mean(axis=0)
plt.plot(errs)
plt.show()
for c in range(ncols):
    ctrs[:,c,1] = errs


envi.save_image('../data/EMIT_WavelengthCenters_20220117.hdr',np.array(ctrs,dtype=np.float32),ext='',force=True)

wvl = ctrs.mean(axis=1)
fwhm = np.zeros(nrows)
chn = np.arange(328)

np.savetxt('../data/EMIT_Wavelengths_20220117.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
