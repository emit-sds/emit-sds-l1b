# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from scipy.linalg import norm, eigh
import sys, os
from r_pca import R_pca
from emit_fpa import linearity_nbasis

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
    parser.add_argument('output')
    args = parser.parse_args()

    data = None
    use = [90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]
    use = [165,240,315,390,840,915,990,1065,1140]

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

    # ignore extrema
    #data[:,42000:]=np.tile(data[:,42000:42001],(1,len(grid[42000:])))
    #data[:,:500]=np.tile(data[:,500:501],(1,500))
    print(data.shape[0],'datapoints')

    for i in range(data.shape[0]):
      data[i,:] = medfilt(data[i,:])

    data[np.logical_not(np.isfinite(data))]=0

    data_resamp = []
    grid_resamp = np.arange(0,2**16,10,dtype=int)
    for d in data:
        data_resamp.append(interp1d(grid,d)(grid_resamp))
    data_resamp = np.array(data_resamp)
   
    model = PCA(n_components=linearity_nbasis)
    model.fit(data_resamp)
    mu_resamp = model.mean_.copy()
    evec_resamp = model.components_.copy()

   #mu = data_resamp.mean(axis=0)
   #C = np.cov(data_resamp-mu,rowvar=False)

   #ev, evec_resamp = eigh(C)
   #evec_resamp = np.fliplr(evec_resamp[:,-linearity_nbasis:]).T
    
    # Resample and renormalize eigenvectors
    evec = []
    for ev_resamp in evec_resamp:
       ev = interp1d(grid_resamp, ev_resamp, bounds_error=False, 
                     fill_value='extrapolate')(grid)
       ev_normalized = ev / norm(ev)
       evec.append(ev_normalized)
    evec = np.array(evec)
    mu = interp1d(grid_resamp, mu_resamp, 
                   bounds_error=False, fill_value='extrapolate')(grid)
    
    combined = np.concatenate((mu[np.newaxis,:],evec),axis=0)
    combined = combined.astype(np.float32)
    envi.save_image(args.output+'.hdr',combined,ext='',force=True)

    plt.plot(combined.T)
    plt.show()

if __name__ == '__main__':

    main()
