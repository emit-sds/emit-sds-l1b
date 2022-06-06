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
c,wl,f = np.loadtxt('../../data/EMIT_Wavelengths_20220421.txt').T 
x_all, y_all, x2_all = [],[],[]

for fieldpoint,color in [(40,[0.2,0.2,0.8]),(340,[0.8,0.2,0.2]),(640,[0.2,0.8,0.2]),(940,[0.8,0.8,0.2]),(1240,[0.8,0.2,0.8])]:
    xc,x,y = [],[],[]
    files = glob.glob('/beegfs/scratch/drt/20211113_EMIT_SRF/all/*Field%i*subframe_darksub_pedestal.txt'%fieldpoint)
    for fil in files:
       print(fil)
       try:
            chn,fwhm = np.loadtxt(fil,skiprows=2).T
            chn = np.array(chn,dtype=int)
            fwhm = np.array(fwhm)
            if chn.size<2:
                continue
            # filter out noisy fits. We know the FWHM is in the 1-1.4 pixel range
            use = np.logical_and(fwhm>1.0,fwhm<1.4)
            chn = chn[use]
            fwhm = fwhm[use]

            plt.plot(wl[chn],fwhm,'.',color=color,alpha=1.0,markersize=5)

            for c,f, in zip(chn,fwhm):
              x.append(c)
              xc.append(wl[c])
              y.append(f)

       except ValueError:
           continue
    x,xc,y = np.array(xc),np.array(x),np.array(y)
    i = np.argsort(x)
    x,y,xc = x[i],y[i],xc[i]
    x_all = np.concatenate([x_all, xc])
    x2_all = np.concatenate([x2_all, np.ones(len(x))*fieldpoint])
    y_all = np.concatenate([y_all, y])

with open('../data/plots/EMIT_SRF_fits.txt','w') as fout:           
   for x,x2,y in zip(x_all,x2_all,y_all):
      print(x,x2,y)
      fout.write('%8.6f %8.6f %8.6f\n'%(x,x2,y))

p  = bisplrep(x_all,x2_all,y_all,kx=1,ky=1,s=10)
znew = bisplev(np.arange(328),np.arange(1280),p)

envi.save_image('../data/../EMIT_SRF_20220422.hdr',np.array(znew,dtype=np.float32),ext='',force=True)

channel_widths = abs(np.diff(wl))
channel_widths = np.concatenate((channel_widths,
                                np.array([channel_widths[-1]])),axis=0)
fwhm = np.mean(znew, axis=1) * channel_widths
np.savetxt('../../data/EMIT_Wavelengths_20220422.txt',
          np.c_[np.arange(len(wl)), wl, fwhm],fmt='%10.8f')

plt.ylim([0,20])   
plt.xlabel('Wavelength (nm)')
plt.ylabel('FWHM')
plt.box(False)
plt.grid(True)

dwl = (wl[0]-wl[1])*1000.0
fig=plt.figure()
ax = plt.axes()
X, Y = np.meshgrid(np.arange(1280),np.arange(328))
im = ax.contour(X,Y,znew*(dwl),inline=True,colors=['k'])
plt.clabel(im, inline=1, fontsize=12, inline_spacing=0, fmt='%1.1f nm ')
plt.xlabel('Cross track position')
plt.ylabel('Spectral channel')
plt.xlim([-50,1290])
plt.show()
