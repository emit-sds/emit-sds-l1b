# David R Thompson
import numpy as np
import os, sys, glob, os.path
import pylab as plt
import numpy as np
from scipy.interpolate import bisplrep,bisplev
from spectral.io import envi
sys.path.append('../utils')
from lowess import lowess
from emit_fpa import first_valid_row, last_valid_row 

# This script combines spectral response function data from all the
# monochromator sweeps.  It combines them into an array and fits 
# smooth functions to this relationship.

directory = os.path.split(os.path.abspath(__file__))[0]
c,wl,f = np.loadtxt('../../data/aviris3/AVIRIS3_Wavelengths_20230609.txt').T 

color = 'k'
files = glob.glob('/beegfs/scratch/drt/20230608_AVIRIS3_SRF/*darksub_pedestal.txt')
x, xwl, y = [],[],[]
for fil in files:
   print(fil)
   chn,fwhm = np.loadtxt(fil,skiprows=2).T
   chn = np.array(chn,dtype=int)
   fwhm = np.array(fwhm)
   if chn.size<2:
       continue

   # filter out noisy fits. We know the FWHM is in the 1-1.15 pixel range
   use = np.logical_and(fwhm>1.0,fwhm<1.15)
   chn = chn[use]
   fwhm = fwhm[use]

   plt.plot(wl[chn],fwhm,'.',color=color,alpha=1.0,markersize=5)

   for c,f, in zip(chn,fwhm):
     x.append(c)
     xwl.append(wl[c])
     y.append(f)


x,xwl,y = np.array(x),np.array(xwl),np.array(y)
i = np.argsort(x)
x,y,xwl = x[i],y[i],xwl[i]

xall = np.arange(len(wl))
p = np.polyfit(x,y,3)
ynew = np.polyval(p,xall)

fits_file = '../../data/aviris3/plots/AVIRIS3_SRF_fits.txt'
np.savetxt(fits_file, np.c_[x,y], fmt='%10.8f')
plt.plot(xwl,y,'k.')
plt.plot(wl,ynew,'r')
plt.ylim([0,2])   
plt.xlabel('Wavelength (nm)')
plt.ylabel('FWHM (pixels)')
plt.box(False)
plt.grid(True)
plt.show()

channel_widths = abs(np.diff(wl))
channel_widths = np.concatenate((channel_widths,
                                np.array([channel_widths[-1]])),axis=0)
fwhm = ynew * channel_widths
np.savetxt('../../data/aviris3/AVIRIS3_Wavelengths_20230610.txt',
          np.c_[xall, wl, fwhm],fmt='%10.8f')

