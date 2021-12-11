# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from copy import deepcopy
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from scipy.optimize import minimize
from scipy.interpolate import BSpline,interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
import json
from scipy.optimize import minimize
from numba import jit
from fixghost import fix_ghost
import ray
import ray.services

rayargs={'num_cpus':40}
ray.init(**rayargs)

eps = 1e-6


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def serialize_ghost_config(config):

  x = [config['blur_spatial'],config['blur_spectral']]
  for i in range(len(config['orders'])):
      x.append(config['orders'][i]['intensity'])
  return x    


def deserialize_ghost_config(x, config):
  ghost_config = deepcopy(config) 
  if (len(x)-2) != len(config['orders']):
    raise IndexError('bad state vector size')
  ghost_config['blur_spatial'] = x[0]
  ghost_config['blur_spectral'] = x[1]
  for i in range(len(config['orders'])):
    ghost_config['orders'][i]['intensity'] = x[2+i]
  return ghost_config   


@ray.remote
def frame_error(frame, new_config):
    fixed = fix_ghost(frame, new_config)
    half = 640
    max_left = np.percentile(frame[:,:half],99)
    max_right = np.percentile(frame[:,half:],99)
    if max_left>max_right:
        return np.mean(pow(fixed[:,half:],2))
    else:
        return np.mean(pow(fixed[:,:half],2))


def err(x, frames, config):
    new_config = deserialize_ghost_config(x, config)
    jobs = [frame_error.remote(frame, new_config) for frame in frames]
    errs = ray.get(jobs)
    print(x,errs)
    return sum(errs)
  

def main():

    description = "Optimize ghost model"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('input',nargs='+')
    parser.add_argument('output')
    args = parser.parse_args()

    frames = []
    for infile in args.input:
        print(infile)
        I = envi.open(infile)
        bip = np.squeeze(I.load())
        bil = bip.T
        frames.append(bil)
    frames = np.array(frames)
    print(frames.shape)
    with open(args.config,'r') as fin:
        ghost_config = json.load(fin)

    x0 = serialize_ghost_config(ghost_config)
    opts = {'max_iters':5}
    best = minimize(err, x0, args=(frames, ghost_config), 
        options=opts, bounds=[(1,10),(1,10)]+[(0,0.01) for q in x0[2:]])
    best_config = deserialize_ghost_config(best.x, ghost_config)
    
    with open(args.output,'w') as fout:
        fout.write(json.dumps(best_config,indent=2))

if __name__ == '__main__':
    main()
