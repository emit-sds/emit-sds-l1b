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


basedir = '/beegfs/scratch/drt/20220112_CWIS2/20211210_lasers/'
filepath = basedir+'20211210_LaserSphere_clip_darksub_pedestal'
I = envi.open(filepath+'.hdr')
lines = int(I.metadata['lines'])

# Laser info 
# Wavelengths from Sven.  We start with an initial first guess of the channels
# based on her mesurements.
wavelengths = np.array([2064.35,1550.6,1064,632.83,532,406.7]) 
uncert = np.array([0.01,0.1,1,0.01,1,0.1]) 
channels = np.array([84,153,218,276,290,307])
nlasers = len(channels)

# Change to refractive wavelength of vacuum
#index_of_refraction = np.array([1.000268,1.000269,1.000271,1.000273,1.000277])
#wavelengths = wavelengths * index_of_refraction
 
# These channels are reported in whole-FPA format, but the EMIT FPA begins reading at row 6
nrows, ncols = 328, 1280

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
       totals[:,use] = frame[:,use]
       counts[:,use] = counts[:,use] + 1.0

frame = totals / counts

margin = 5
observed = [[] for c in channels]
for col in range(20,ncols-20):
    for i, w, chn in zip(range(nlasers), wavelengths, channels):
            idx = np.arange(int(chn-margin),int(chn+margin+1), dtype=int)
            row,_,_ = find_peak(frame[idx,int(round(col))])
            row = row+idx[0] 
            observed[i].append([col,row,w])

x = np.zeros((nlasers,ncols))
y = np.zeros((nlasers,ncols))
for i in range(len(observed)):
    D = np.array(observed[i])
    use = np.logical_and(D[:,0]>200,D[:,0]<1000)
    p = np.polyfit(D[use,0],D[use,1],1)
    resid = D[:,1]-np.polyval(p, D[:,0])
    #use = resid < np.median(resid)
    #p = np.polyfit(D[use,0],D[use,1],1)
    x[i,:] = np.polyval(p, np.arange(ncols))
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


envi.save_image('../data/CWIS_WavelengthCenters_20220331.hdr',np.array(ctrs,dtype=np.float32),ext='',force=True)

