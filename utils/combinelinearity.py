# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
import sys, os


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
    parser.add_argument('--nev',default=5)
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
    data = np.array(data) #[/ np.arange(data.shape[1])
    mu = data.mean(axis=0)
    zm = data - mu
    use = np.arange(100,40000,10)
    plt.plot(np.std(zm[:,use],axis=0))
    
    use = np.arange(0,47000,10)
    C = np.cov(zm[:,use], rowvar=False)
    C[np.logical_not(np.isfinite(C))]=0
    ev,vec = np.linalg.eig(C)
    resamp = []
    for v in vec[:,:args.nev].T:
       resamp.append(interp1d(use,v,fill_value='extrapolate',bounds_error=False)(np.arange(2**16)))
    resamp = np.array(resamp)
    envi.save_image(args.output+'.hdr',np.concatenate((mu[np.newaxis,:],resamp),axis=0),ext='',force=True)

if __name__ == '__main__':

    main()
