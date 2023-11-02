# David R Thompson
import numpy as np
import pylab as plt
from scipy.interpolate import splev, splrep

# Load wavelengths and fwhm data
c, wl, _ = np.loadtxt('../../data/cpm/CPM_Wavelengths_20231001.txt').T *1000
c, wlb, field30, field175, field320, field465, field610_= \
           np.loadtxt('../../data/cpm/ancillary/CPM_FWHM_Measurements.txt').T

# Choose the center field point
y = field320

# Fit a smoothing spline
valid = np.where(np.logical_and(y>1.0,y<1.25))[0]
tck = splrep(c[valid],y[valid],s=0.5)
smoothed = splev(c[valid],tck)

# Remove monochromator outliers and refit
inliers = valid[abs(y[valid]-smoothed)<(np.std(y[valid]-smoothed))]
tck = splrep(c[inliers],y[inliers],s=1)
filtered = splev(c,tck)

if False:
    plt.plot(wl[valid],y[valid],'r.')
    plt.plot(wl[inliers],y[inliers],'k.')
    plt.plot(wl[valid],smoothed,'r')
    plt.plot(wl,filtered,'k')
    plt.ylim([0,2])
    plt.box(False)
    plt.grid(True)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('FWHM (channels)')
    plt.show()

# Translate channels to nm
del_wl = np.abs(np.diff(wl))
del_wl = np.r_[del_wl, del_wl[-1]]
fwhm = filtered * del_wl

# Write out, converting from nm to microns 
np.savetxt('../../data/cpm/CPM_Wavelengths_20231015.txt', 
    np.c_[c,wl/1000.0,fwhm/1000.0], fmt='%10.8f')
