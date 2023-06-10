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
from scipy.ndimage import binary_dilation
from numpy.random import randn
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import make_pipeline
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


# We will analyze the laser sphere data.  It has already
# been clipped to the standard subframe, dark-subtracted, and 
# pedestal-shift corrected
basedir = '/beegfs/scratch/drt/20230608_AVIRIS3_LaserSphere/'
filepath = basedir+'20221212_Laser_Sphere_darksub_pedestal'

I = envi.open(filepath+'.hdr')
lines = int(I.metadata['lines'])

#We start with an initial first guess of the channels
wavelengths = np.array([2064.35,1550.6,1064,632.83,532,406.7]) 
uncert = np.array([0.01,0.1,1,0.01,1,0.1]) 
channels = np.array([83,152,218,276,289,306])
nlasers = len(channels)

# Change to refractive wavelength of vacuum - already done?
 
# These channels are reported in whole-FPA format, but the EMIT FPA begins reading at row 7
nrows, ncols = 328, 1280
margin = 5

totals = np.zeros((nrows,ncols),dtype=np.float32)
counts = np.zeros((nrows,ncols),dtype=np.float32)
with open(filepath,'rb') as fin:
   for line in range(lines):
       frame = np.fromfile(fin,dtype=np.float32,count=ncols*nrows)
       frame = frame.reshape((nrows,ncols))
       magnitude = frame.mean(axis=0) 
      #plt.plot(magnitude)
      #plt.show()
       use = magnitude>48
       totals[:,use] = totals[:,use] + frame[:,use]
       counts[:,use] = counts[:,use] + 1.0
frame = totals / counts

# our list of laser fits, one sublist per laser
observed = [[] for c in channels]

for col in range(40,ncols-40):
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
    np.savetxt('../../data/aviris3/plots/AVIRIS3_Laser_%i_ColRow.txt'%wavelengths[i],D,fmt='%8.6f')
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
for i in [640]:
  print('column',i)
  x0 = np.array([2645,1])
  best = minimize(errs,x0,args=(y[:,i],x[:,i])) 
  p = np.polyfit(x[:,i],y[:,i],1)
  ctrs[:,i,0] = gen_wavelengths(best.x)

# save out wavelength center matrix
envi.save_image('../../data/aviris3/AVIRIS3_WavelengthCenters_20230609.hdr',np.array(ctrs,dtype=np.float32),ext='',force=True)

# save out ascii text file, no FWHMs, in microns
wvl = np.squeeze(ctrs[:,:,0]).mean(axis=1)
wvl = ctrs[:,640,0]
fwhm = np.zeros(nrows)
chn = np.arange(328)
np.savetxt('../../data/aviris3/AVIRIS3_Wavelengths_20230609.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
