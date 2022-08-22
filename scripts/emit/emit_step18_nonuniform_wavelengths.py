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
from scipy.signal import medfilt
from numpy.random import randn
import json

# Fit a gaussian to a single peak in an ordered series of data
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


# We will analyze the laser sphere data from TVAC2.  It has already
# been clipped to the standard subframe, dark-subtracted, and 
# pedestal-shift corrected
basedir = '/beegfs/scratch/drt/20220124_EMIT_LaserSphere/'
I = envi.open(basedir+'20211112_211905_UTC_LaS_Fields-40-1319_subframe_darksub_pedestal.hdr')
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
 
# These channels are reported in whole-FPA format, but the EMIT FPA begins reading at row 7
row_offset = 7
nrows, ncols = 328, 1280
channels = channels - row_offset 

margin = 5

# our list of laser fits, one sublist per laser
observed = [[] for c in channels]

# Find a the spatial location of the lasers, and fit each laser peak
# put them as column, row, wavelength triplets into the "observed" list
for line in range(I.shape[0]):
   frame = np.squeeze(I[line,:,:]).T 
   col,amp,_ = find_peak(medfilt(frame.mean(axis=0),5))
   print(col,amp)
   if amp<100:
       continue
   for i, w, chn in zip(range(nlasers), wavelengths, channels):
       idx = np.arange(int(chn-margin),int(chn+margin+1), dtype=int)
       row,_,_ = find_peak(frame[idx,int(round(col))])
       row = row+idx[0] 
       observed[i].append([col,row,w])

# Now plot the result, and save the laser fits to a file for later plotting
x = np.zeros((nlasers,ncols))
y = np.zeros((nlasers,ncols))
for i in range(len(observed)):
    D = np.array(observed[i])
    # map row to column for a single laser
    p = np.polyfit(D[:,0],D[:,1],1)
    x[i,:] = np.polyval(p,np.arange(ncols))
    y[i,:] = wavelengths[i]
    plt.plot(D[:,0],D[:,1]-p[-1],'.')
    plt.plot(np.arange(ncols),x[i,:]-p[-1],'k')
    np.savetxt('../../data/plots/EMIT_Laser_%i_ColRow.txt'%wavelengths[i],D,fmt='%8.6f')
plt.show()

# now we fit the nonuniform dispersion curve
model_wavelengths=[380.00000,439.00000,498.00000,557.00000,616.00000,675.00000,800.00000,925.00000,1050.00000,1175.00000,1300.00000,1542.00000,1784.00000,2026.00000,2268.00000,2510.00000]
model_dispersions=[7.39600,7.41900,7.43300,7.44200,7.44800,7.45200,7.45700,7.45800,7.45600,7.45300,7.45000,7.43900,7.42600,7.41100,7.39300,7.37300]
average = 7.431625

# FPA is flipped
model_dispersions = -np.array(model_dispersions)

dispersion = interp1d(model_wavelengths, model_dispersions, bounds_error=False, fill_value='extrapolate')

def gen_wavelengths(x, model_wavelengths=model_wavelengths, model_dispersions=model_dispersions):
  start = x[0]
  stretch = x[1]
  wls = [start]
  for i in range(1,nrows):
     wls.append(wls[-1]+dispersion(wls[-1])*stretch)
  return np.array(wls)

def errs(x,actual_wavelengths,actual_columns):
  wl = gen_wavelengths(x)
  predicted_wl = interp1d(np.arange(len(wl)), wl, bounds_error=False, fill_value='extrapolate')(actual_columns)
  err = np.sum((actual_wavelengths - predicted_wl)**2)
  return err 


ctrs = np.zeros((nrows, ncols, 2))

# Perform the fit for each column
ctrs = np.zeros((nrows, ncols, 2))
for i in [640]:#range(ncols):
  print('column',i)
  x0 = np.array([2645,1])
  best = minimize(errs,x0,args=(y[:,i],x[:,i])) 
  p = np.polyfit(x[:,i],y[:,i],1)
  ctrs[:,i,0] = gen_wavelengths(best.x)

# save out wavelength center matrix
envi.save_image('../../data/EMIT_WavelengthCenters_20220805.hdr',np.array(ctrs,dtype=np.float32),ext='',force=True)

# save out ascii text file, no FWHMs, in microns
wvl = np.squeeze(ctrs[:,:,0]).mean(axis=1)
wvl = ctrs[:,640,0]
fwhm = np.zeros(nrows)
chn = np.arange(328)
np.savetxt('../../data/EMIT_Wavelengths_20220805.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
