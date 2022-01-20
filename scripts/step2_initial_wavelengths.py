# David R Thompson
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d

# Laser info from 20211112_211905_UTC_LaS_Fields-40-1319
# Wavelengths from Christine
wavelengths = np.array([1950,1064,633,532,405]) 
channels = np.array([100.36,219.29,277.15,290.62,307.66])

# These channels are reported in whole-FPA format, but the EMIT FPA begins reading at row 6
row_offset = 6
nrows = 328
channels = channels - row_offset 

# Change to refractive wavelength of vacuum
index_of_refraction = np.array([1.000268,1.000269,1.000271,1.000273,1.000277])
wavelengths = wavelengths * index_of_refraction
 
# Calculate the polynomial fit to all lasers
chn = np.arange(328)
p = np.polyfit(channels,wavelengths,1)
wvl = np.polyval(p, chn)

fwhm = np.zeros(328)

np.savetxt('../data/EMIT_Wavelengths_20220117.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
