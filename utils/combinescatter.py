# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from scipy.stats import norm
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from scipy.linalg import inv
import sys, os
from lowess import lowess

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def sum_of_gaussians(x, mean1=0, amplitude1=1., sigma1=1.,
                        amplitude2=1., sigma2=1.,
                        amplitude3=1., sigma3=1.):
    return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) +
            amplitude2 * np.exp(-0.5 * ((x - mean1) / sigma2)**2)  +
            amplitude3 * np.exp(-0.5 * ((x - mean1) / sigma3)**2))


left, right, shortw, longw = 25, 1265, 21, 314


def main():

    description = "Calculate linearity basis"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--spatial',action='store_true')
    parser.add_argument('--manual',default=1.0)
    parser.add_argument('output')
    args = parser.parse_args()

    ctr, mean1, amp1, sigma1, amp2, sigma2, amp3, sigma3, err = np.loadtxt(args.input).T

    if False:
        plt.plot(ctr,sigma1,'k.')
        plt.plot(ctr,amp1,'r.')
        
        plt.figure()
        plt.plot(ctr,sigma2,'k.')
        plt.plot(ctr,amp2,'r.')
        
        plt.figure()
        plt.plot(ctr,sigma3,'k.')
        plt.plot(ctr,amp3,'r.')
        
        plt.show()

     

    if args.spatial:

        use = np.logical_and(amp2>0,amp3>0)
        sigma1 = np.median(sigma1[use])
        amp2=np.median(amp2[use])
        sigma2 = np.median(sigma2[use])
        amp3 = np.median(amp3[use])
        sigma3 = np.median(sigma3[use])
        L = np.zeros((1280,1280))
        x = np.arange(1280)

        # set amplitude relative to main peak with area of unity
        ampn = norm.pdf(x,640,sigma1)
        ampn = ampn / ampn.sum()
        print('amplitude',ampn.max())
        amp2 = amp2 * ampn.max() * float(args.manual)
        amp3 = amp3 * ampn.max() * float(args.manual)
        print('spatial - using',amp2,sigma2,amp3,sigma3)

        for i in range(left, right+1):
             if i%100==0:
                 print(i,'/',1280)
             # sigma of 0.5 is arbitrary because central peak magnitude is zero
             L[i,:] = sum_of_gaussians(x,i,0,0.5,amp2,sigma2,amp3,sigma3)
             L[i,:left] = 0
             L[i,(right+1):] = 0
             
        L = L+np.eye(1280)

        # Normalize each row to conserve photons
        for i in range(left, right+1):
             L[i,:] = L[i,:] / np.sum(L[i,:])

    else:

        use = np.logical_and(amp2>0,amp3>0)
        sigma1 = np.median(sigma1[use])
        amp2=np.median(amp2[use])
        sigma2 = np.median(sigma2[use])
        amp3 = np.median(amp3[use])
        sigma3 = np.median(sigma3[use])
        L = np.zeros((480,480))
        x = np.arange(480)
  
        # set amplitude relative to main peak with area of unity
        print(sigma1.shape)
        ampn = norm.pdf(x,240,sigma1)
        ampn = ampn / ampn.sum()
        print('amplitude',ampn.max())
        amp2 = amp2 * ampn.max() * float(args.manual)
        amp3 = amp3 * ampn.max() * float(args.manual)
        print('spectral - using',amp2,sigma2,amp3,sigma3)

        for i in range(shortw,longw+1):
             if i%100==0:
                 print(i,'/',480)
             # sigma of 0.5 is arbitrary because central peak magnitude is zero
             L[i,:] = sum_of_gaussians(x,i,0,0.5,amp2,sigma2,amp3,sigma3)
             L[i,:shortw] = 0
             L[i,(longw+1):] = 0

        L = L+np.eye(480)

        # Normalize each row to conserve photons
        for i in range(shortw, longw+1):
             L[i,:] = L[i,:] / np.sum(L[i,:])
  
    Linv = inv(L)
    Linv = Linv.astype(np.float32)
    envi.save_image(args.output+'.hdr',Linv,ext='',force=True)



if __name__ == '__main__':

    main()
