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

eps = 1e-4
optimize_blur=False


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def serialize_ghost_config(config):

  if optimize_blur:
      x = [config['blur_spatial'],config['blur_spectral']]
  else:
      x = []
  for i in range(len(config['orders'])):
      x.append(np.log(config['orders'][i]['intensity']))
  return x    


def deserialize_ghost_config(x, config):
  ghost_config = deepcopy(config) 
  if optimize_blur:
    if (len(x)-2) != len(config['orders']):
      raise IndexError('bad state vector size')
    ghost_config['blur_spatial'] = x[0]
    ghost_config['blur_spectral'] = x[1]
    ind = 2
  else:
    if (len(x)) != len(config['orders']):
      raise IndexError('bad state vector size')
    ind = 0
  for i in range(len(config['orders'])):
    ghost_config['orders'][i]['intensity'] = np.exp(x[ind+i])
  return ghost_config   


#@ray.remote
def frame_error(frame, new_config):
    fixed = fix_ghost(frame, new_config)
    half = 640
    max_left = np.percentile(frame[:,:half],99)
    max_right = np.percentile(frame[:,half:],99)
    if max_left>max_right:
        return np.mean(pow(fixed[:,half:],2)) / np.mean(pow(frame[:,half:],2)) 
    else:
        return np.mean(pow(fixed[:,:half],2)) / np.mean(pow(frame[:,:half],2))


def err(x, frames, ghost_config):
    new_config = deserialize_ghost_config(x, ghost_config)
    jobs = [frame_error(frame, new_config) for frame in frames]
    errs = np.array(jobs)#ray.get(jobs)
    disp = x.copy()
    if optimize_blur:
        disp[2:] = np.exp(x[2:])
    else:
        disp = np.exp(x)
    print(disp,errs)
    return sum(errs)
 
@ray.remote
def partial(x, i, frames, ghost_config, base_cost):
    x_perturb = x.copy()
    x_perturb[i] = x[i] + eps
    perturb_cost = err(x_perturb, frames, ghost_config)
    return (perturb_cost - base_cost)/eps


def jac(x, frames, ghost_config):
    base_cost = err(x,frames,ghost_config)
    jobs  = [partial.remote(x,i,frames,ghost_config,base_cost) for i in range(len(x))]
    derivs = ray.get(jobs)
    return np.array(derivs)


def main():

    description = "Optimize ghost model"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('input',nargs='+')
    parser.add_argument('output')
    args = parser.parse_args()

    frames = []
    for infile in args.input:
        I = envi.open(find_header(infile))
        bip = np.squeeze(I.load())
        bil = bip.T
        frames.append(bil)
    frames = np.array(frames)
    with open(args.config,'r') as fin:
        ghost_config = json.load(fin)
 
    x0 = serialize_ghost_config(ghost_config)
    opts = {'max_iters':10}
    best = minimize(err, x0, args=(frames, ghost_config), jac=jac,
        options=opts)#, bounds=[(1,10),(1,10)]+[(0,0.01) for q in x0[2:]])
    best_config = deserialize_ghost_config(best.x, ghost_config)
    
    with open(args.output,'w') as fout:
        fout.write(json.dumps(best_config,indent=2))

if __name__ == '__main__':
    main()
