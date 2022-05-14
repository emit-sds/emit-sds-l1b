# David R Thompson
import numpy as np
import pylab as plt
import os, sys
from spectral.io import envi
from scipy.interpolate import interp1d
from PIL import Image

c,wl,fwhm = np.loadtxt('../../data/emit/EMIT_Wavelengths_20220422.txt').T * 1000.0
c,rcc,unc = np.loadtxt('../../data/emit/EMIT_RadiometricCoeffs_20220504.txt').T
swl, irr = np.loadtxt('../../data/emit/kurudz_0.1nm.dat').T
irr = irr * 1e3 / 1e4 # convert mW / (m2 nm) to uW / (cm2 nm)
irr = interp1d(swl,irr,bounds_error=False,fill_value='extrapolate')(wl)

if False:
  im = envi.open("/beegfs/scratch/drt/20220513_SaturationMap/Landsat8_Global_Median_2021_resampled.hdr")
  lwl = np.array([440,480,565,665,865,1610,2200])
  lb = np.array([np.argmin(abs(b-wl)) for b in lwl])
  im = im.load()
   
  # Calculate solar zenith associated with saturation radiance
  dn = 48000
  zeniths = -np.ones(im.shape)
  for r in range(im.shape[0]):
    print(r)
    for c in range(im.shape[1]):
      for b in range(im.shape[2]):
        satlevel = dn / (im[r,c,b]  * irr[lb[b]] / np.pi / rcc[lb[b]])
        if satlevel>1: 
           continue
        zen = np.arccos(satlevel)
        zen = zen / 2.0 / np.pi * 360.0
        zeniths[r,c,b] = zen
  zeniths = np.nanmax(zeniths,axis=2)

  envi.save_image('"/beegfs/scratch/drt/20220513_SaturationMap/max_zenith.hdr',zeniths,ext='',force=True)
      

# Now calculate cloud screening parameters
dark = 8300 # actually ranges from 7800-8300 DN
d1,d2,d3 = dark,dark,dark
c1 = np.argmin(abs(wl-450))
c2 = np.argmin(abs(wl-1250))
c3 = np.argmin(abs(wl-1650))

alpha1 = 0.51
alpha2 = 0.56
alpha3 = 0.29
b1 = irr[c1] * alpha1 / np.pi / rcc[c1]
b2 = irr[c2] * alpha2 / np.pi / rcc[c2]
b3 = irr[c3] * alpha3 / np.pi / rcc[c3]

print('c1',c1)
print('c2',c2)
print('c3',c3)
print('b1',b1)
print('b2',b2)
print('b3',b3)
print('d1',d1)
print('d2',d2)
print('d3',d3)
