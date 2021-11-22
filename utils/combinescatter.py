# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from scipy.linalg import norm,inv
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
            amplitude2 * np.exp(-0.5 * ((x - mean1) / sigma2)**2) +
            amplitude3 * np.exp(-0.5 * ((x - mean1) / sigma3)**2))




def main():

    description = "Calculate linearity basis"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--spatial',action='store_true')
    parser.add_argument('output')
    args = parser.parse_args()

    ctr, mean1, amp1, sigma1, amp2, sigma2, amp3, sigma3, err = np.loadtxt(args.input).T

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

        amp2=np.median(amp2)
        sigma2 = np.median(sigma2)
        amp3 = np.median(amp3)
        sigma3 = np.median(sigma3)
        L = np.zeros((1280,1280))
        x = np.arange(1280)
        for i in range(1280):
             if i%100==0:
                 print(i,'/',1280)
             L[i,:] = sum_of_gaussians(x,i,0,1,amp2,sigma2,amp3,sigma3)

    else:

        use = np.logical_and(amp2>0,amp3>0)
        amp2=np.median(amp2[use])
        sigma2 = np.median(sigma2[use])
        amp3 = np.median(amp3[use])
        sigma3 = np.median(sigma3[use])
        L = np.zeros((480,480))
        x = np.arange(480)
        for i in range(480):
             if i%100==0:
                 print(i,'/',480)
             L[i,:] = sum_of_gaussians(x,i,0,1,amp2,sigma2,amp3,sigma3)
  
    Linv = inv(L)
    Linv = Linv.astype(np.float32)
    envi.save_image(args.output+'.hdr',Linv,ext='',force=True)



    
  # for fi,infilepath in enumerate(args.input):
  #     print(fi,'/',len(args.input))
  #     if not any([str(u) in infilepath for u in use]):
  #         continue

  #     x = envi.open(infilepath+'.hdr').load()      
  #     x = np.squeeze(x)
  #     if data is None:
  #         data = x
  #     else:
  #         data = np.concatenate((data,x),axis=0)

  # grid = np.arange(2**16, dtype=float)
  # data = np.array(data) 
  # data[0] = 0
  # data[np.logical_not(np.isfinite(data))]=0

  # # remove large outlier values
  # norms = norm(data,axis=1)
  # mu = np.mean(norms)
  # std = np.std(norms)
  # use = (norms-mu)<(std*3)
  # data = data[use,:]
  # print(data.shape[0],'datapoints')

  # for i in range(data.shape[0]):
  #   data[i,:] = medfilt(data[i,:])

  # data[np.logical_not(np.isfinite(data))]=0

  # pca = PCA(args.nev)
  # pca.fit(data)
  # resamp = pca.components_
  # mu = pca.mean_
  # 
  ##mu = data.mean(axis=0)
  ##plt.plot(mu)
  ##plt.show()
  ##zm = data - mu
  ##
  ##use = np.arange(0,47000,10)
  ##C = np.cov(zm[:,use], rowvar=False)
  ##C[np.logical_not(np.isfinite(C))]=0

  ##print('Calculating eigenvectors')
  ##ev,vec = np.linalg.eigh(C)
  ##
  ##print(ev[-args.nev:])
  ##print(mu)
  ##resamp = []
  ##for i in range(args.nev):
  ##   v = vec[:,-(i+1)]
  ##   #plt.plot(v)
  ##   v = interp1d(use,v,fill_value='extrapolate',bounds_error=False)(np.arange(2**16))
  ##   v = v / norm(v)
  ##   resamp.append(v)
  ###plt.show()

  ##resamp = np.array(resamp)


if __name__ == '__main__':

    main()
