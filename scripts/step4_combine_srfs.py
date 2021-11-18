# David R Thompson
import numpy as np
import os, sys, glob
import pylab as plt
import numpy as np
from scipy.interpolate import bisplrep,bisplev
from spectral.io import envi
sys.path.append('../utils')
from lowess import lowess


c,wl,f = np.loadtxt('/home/drt/src/emit-sds-l1b/data/EMIT_Wavelengths_20211104.txt').T
x_all, y_all, x2_all = [],[],[]
for fieldpoint,color in [(40,'r'),(340,'b'),(640,'g'),(940,'y'),(1240,'c')]:
    xc,x,y = [],[],[]
    files = glob.glob('/beegfs/scratch/drt/20211113_EMIT_SRF/all/*Field%i*'%fieldpoint)
    for fil in files:
        print(fil)
        try:
            chn,fwhm = np.loadtxt(fil).T
            chn = np.array(chn,dtype=int)
            fwhm = np.array(fwhm)
            if chn.size<2:
                continue
            plt.plot(wl[chn]*1000.0,fwhm,'.',color=color,alpha=0.3,markersize=10)
            for c,f, in zip(chn,fwhm):
              x.append(c)
              xc.append(wl[c])
              y.append(f)
        except ValueError:
            continue
    x,xc,y = np.array(xc),np.array(x),np.array(y)
    i = np.argsort(x)
    x,y,xc = x[i],y[i],xc[i]
    #ys = lowess(xc,y,f=0.1)
    x_all = np.concatenate([x_all, xc])
    x2_all = np.concatenate([x2_all, np.ones(len(x))*fieldpoint])
    #y_all = np.concatenate([y_all, ys])
    y_all = np.concatenate([y_all, y])
    #plt.plot(x,ys,color=color)
           
for x,x2,y in zip(x_all,x2_all,y_all):
  print(x,x2,y)

p  = bisplrep(x_all,x2_all,y_all,kx=1,ky=1)
znew = bisplev(np.arange(480),np.arange(1280),p)

envi.save_image('../data/EMIT_SRF_20211114.hdr',np.array(znew,dtype=np.float32),ext='',force=True)

channel_widths = abs(np.diff(wl))
channel_widths = np.concatenate((channel_widths,
                                np.array([channel_widths[-1]])),axis=0)
fwhm = np.mean(znew, axis=1) * channel_widths
np.savetxt('../data/EMIT_Wavelengths_20211117.txt',
          np.c_[np.arange(len(wl)), wl, fwhm],fmt='%10.8f')

plt.ylim([0,20])   
plt.xlabel('Wavelength (nm)')
plt.ylabel('FWHM')
plt.box(False)
plt.grid(True)

fig=plt.figure()
ax = plt.axes()
im = ax.imshow(znew)
cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
plt.show()
