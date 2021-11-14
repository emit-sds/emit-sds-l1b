# David R Thompson
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d

# laser info
wavelengths = np.array([1956.1,1071.4,640.2,536.1,409.8]) 
channels = np.array([100.36,219.29,277.15,290.62,307.66])
 
chn = np.arange(480)
p = interp1d(channels,wavelengths,fill_value='extrapolate',bounds_error=False)
wvl = p(chn)
fwhm = np.zeros(480)

np.savetxt('../data/EMIT_Wavelengths_20211104.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
