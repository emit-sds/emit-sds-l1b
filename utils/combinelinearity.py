# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from sklearn.decomposition import PCA
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
    parser.add_argument('--nev',type=int,default=2)
    parser.add_argument('output')
    args = parser.parse_args()

    data = None
    use = [90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]

    for fi,infilepath in enumerate(args.input):
        print(fi,'/',len(args.input))
        if not any([str(u) in infilepath for u in use]):
            continue

        x = envi.open(infilepath+'.hdr').load()      
        x = np.squeeze(x)
        if data is None:
            data = x
        else:
            data = np.concatenate((data,x),axis=0)

    grid = np.arange(2**16, dtype=float)
    data = np.array(data) 
    data[0] = 0
    data[np.logical_not(np.isfinite(data))]=0

    # remove large outlier values
    norms = norm(data,axis=1)
    mu = np.mean(norms)
    std = np.std(norms)
    use = (norms-mu)<(std*3)
    data = data[use,:]
    print(data.shape[0],'datapoints')

    for i in range(data.shape[0]):
      data[i,:] = medfilt(data[i,:])

    data[np.logical_not(np.isfinite(data))]=0

    pca = PCA(args.nev)
    pca.fit(data)
    resamp = pca.components_
    mu = pca.mean_
    
    combined = np.concatenate((mu[np.newaxis,:],resamp),axis=0)
    combined = combined.astype(np.float32)
    envi.save_image(args.output+'.hdr',combined,ext='',force=True)

if __name__ == '__main__':

    main()
