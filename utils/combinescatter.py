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
from sklearn import linear_model

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
    i = np.argsort(ctr)
    for field in [mean1, amp1, sigma1, amp2, sigma2, amp3, sigma3, err]:
        field = field[i]
    
    # Hold all spectral points 
    grid = np.arange(480)

    if True:
        use = np.logical_and(amp2>0,amp3>0)
        ransac = linear_model.RANSACRegressor()
        ransac.fit(ctr[use,np.newaxis],sigma1[use])
        sigma1_smoothed = ransac.predict(ctr[:,np.newaxis])
        ransac.fit(ctr[use,np.newaxis],amp1[use])
        amp1_smoothed = ransac.predict(ctr[:,np.newaxis])
        ransac.fit(ctr[use,np.newaxis],sigma2[use])
        sigma2_smoothed = ransac.predict(ctr[:,np.newaxis])
        ransac.fit(ctr[use,np.newaxis],amp2[use])
        amp2_smoothed = ransac.predict(ctr[:,np.newaxis])
        ransac.fit(ctr[use,np.newaxis],sigma3[use])
        sigma3_smoothed = ransac.predict(ctr[:,np.newaxis])
        ransac.fit(ctr[use,np.newaxis],amp3[use])
        amp3_smoothed = ransac.predict(ctr[:,np.newaxis])
    else:
        ctrs, sigma1_smoothed = lowess(sigma1,ctr).T
        ctrs, amp1_smoothed = lowess(amp1,ctr).T
        ctrs, sigma2_smoothed = lowess(sigma2,ctr).T
        ctrs, amp2_smoothed = lowess(amp2,ctr).T
        ctrs, sigma3_smoothed = lowess(sigma3,ctr).T
        ctrs, amp3_smoothed = lowess(amp3,ctr).T

    plt.plot(ctr,sigma1,'k.')
    plt.plot(ctr,amp1,'r.')
    amp1 = interp1d(ctr, amp1_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    sigma1 = interp1d(ctr, sigma1_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    plt.plot(grid, amp1,'r')
    plt.plot(grid, sigma1,'k')
    
    plt.figure()
    plt.plot(ctr,sigma2,'k.')
    plt.plot(ctr,amp2,'r.')
    amp2 = interp1d(ctr, amp2_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    sigma2 = interp1d(ctr, sigma2_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    plt.plot(grid, amp2,'r')
    plt.plot(grid, sigma2,'k')

    plt.figure()
    plt.plot(ctr,sigma3,'k.')
    plt.plot(ctr,amp3,'r.')
    amp3 = interp1d(ctr, amp3_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    sigma3 = interp1d(ctr, sigma3_smoothed, fill_value='extrapolate',bounds_error=False)(grid)
    plt.plot(grid, amp3,'r')
    plt.plot(grid, sigma3,'k')
    plt.show()

    if args.spatial:

        L = np.zeros((1280,1280))
        x = np.arange(1280)

        use = np.logical_and(amp2>0,amp3>0)
        s1 = np.mean(sigma1[use])
        a2 = np.mean(amp2[use])
        s2 = np.mean(sigma2[use])
        a3 = np.mean(amp3[use])
        s3 = np.mean(sigma3[use])

        # set amplitude relative to main peak with area of unity
        ampn = norm.pdf(x,640,s1)
        ampn = ampn / ampn.sum()
        a2 = a2 * ampn.max() * float(args.manual)
        a3 = a3 * ampn.max() * float(args.manual)

        for i in range(left, right+1):
             if i%100==0:
                 print(i,'/',1280)

             # sigma of 0.5 is arbitrary because central peak magnitude is zero
             L[i,:] = sum_of_gaussians(x,i,0,0.5,a2,s2,a3,s3)
             L[i,:left] = 0
             L[i,(right+1):] = 0
             
        L = L+np.eye(1280)

        # Normalize each row to conserve photons
        for i in range(left, right+1):
             L[i,:] = L[i,:] / np.sum(L[i,:])

    else:

        L = np.zeros((480,480))
        x = np.arange(480)
  
        for i in range(shortw,longw+1):

             s1, a2, s2, a3, s3 = sigma1[i],amp2[i],sigma2[i],amp3[i],sigma3[i]
             # set amplitude relative to main peak with area of unity
             ampn = norm.pdf(x,240,s1)
             ampn = ampn / ampn.sum()
             a2 = a2 * ampn.max() * float(args.manual)
             a3 = a3 * ampn.max() * float(args.manual)

             # sigma of 0.5 is arbitrary because central peak magnitude is zero
             L[i,:] = sum_of_gaussians(x,i,0,0.5,a2,s2,a3,s3)
             L[i,:shortw] = 0
             L[i,(longw+1):] = 0
                               
             if i%100==0:
                 print(i,'/',480)
                 plt.semilogy(x,L[i,:])
                 plt.show()

        L = L+np.eye(480)

        # Normalize each row to conserve photons
        for i in range(shortw, longw+1):
             L[i,:] = L[i,:] / np.sum(L[i,:])
  
    Linv = inv(L)
    Linv = Linv.astype(np.float32)
    envi.save_image(args.output+'.hdr',Linv,ext='',force=True)



if __name__ == '__main__':

    main()
