# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from scipy.linalg import norm
import sys, os
from r_pca import R_pca

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')




def main():

    description = "Calculate linearity basis"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('--nev',type=int,default=5)
    parser.add_argument('output')
    args = parser.parse_args()

    data = None
    use = [90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]

    for fi,infilepath in enumerate(args.input):
        if not any([str(u) in infilepath for u in use]):
            continue

        x = envi.open(infilepath+'.hdr').load()      
        x = np.squeeze(x)
        if data is None:
            data = x
        else:
            data = np.concatenate((data,x),axis=0)

    grid = np.arange(2**16, dtype=float)
    data = np.array(data) / grid 
    data[0] = 0
    data[np.logical_not(np.isfinite(data))]=0

    # remove large outlier values
    norms = norm(data,axis=1)
    mu = np.mean(norms)
    std = np.std(norms)
    use = (norms-mu)<(std*3)
    data = data[use,:]

    for i in range(data.shape[0]):
      data[i,:] = medfilt(data[i,:])


    # Robust PCA
    # Implementation of https://arxiv.org/pdf/0912.3599.pdf
    # rpca = R_pca(data.T)
    # data, S = rpca.fit(max_iter=10000, iter_print=100)
    # data = data.T

    data[np.logical_not(np.isfinite(data))]=0
    mu = data.mean(axis=0)
    plt.plot(mu)
    plt.show()
    zm = data - mu
    
    use = np.arange(0,47000,10)
    C = np.cov(zm[:,use], rowvar=False)
    C[np.logical_not(np.isfinite(C))]=0
    ev,vec = np.linalg.eig(C)
    
    ind = np.argsort(ev)
    print(ev[ind[-10:]])
    print(mu)
    resamp = []
    for i in ind[-args.nev:]:
       v = vec[:,i]
      #plt.plot(v)
      #plt.show()
       v = interp1d(use,v,fill_value='extrapolate',bounds_error=False)(np.arange(2**16))
       v = v / norm(v)
       resamp.append(v)

    resamp = np.array(resamp)

    combined = np.concatenate((mu[np.newaxis,:],resamp),axis=0)
    combined = combined.astype(np.float32)
    envi.save_image(args.output+'.hdr',combined,ext='',force=True)

if __name__ == '__main__':

    main()
