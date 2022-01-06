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
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
import json
from scipy.optimize import minimize
from numba import jit
from fixghost import fix_ghost_matrix
from fixghostraster import build_ghost_matrix
import ray
import ray.services

rayargs={'num_cpus':40}
ray.init(**rayargs)


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def serialize_ghost_config(config, coarse):
  x = []
  if not coarse:
      x.append(np.log(config['blur_spectral']))
      x.append(np.log(config['blur_spatial']))
     #for i in range(len(config['orders'])):
     #    x.append(config['orders'][i]['intensity_slope'])
     #    x.append(config['orders'][i]['intensity_offset'])
  else:
      for i in range(len(config['orders'])):
          x.append(np.log(config['orders'][i]['scaling']))
  return x    


def deserialize_ghost_config(x, config, coarse):
  ghost_config = deepcopy(config) 
  ind = 0
  if not coarse:
      ghost_config['blur_spectral'] = np.exp(x[0])
      ghost_config['blur_spatial'] = np.exp(x[1])
     #for i in range(len(config['orders'])):
     #      ghost_config['orders'][i]['intensity_slope'] = x[2+i]
     #      ghost_config['orders'][i]['intensity_offset'] = x[2+i]
  else:
      for i in range(len(config['orders'])):
          ghost_config['orders'][i]['scaling'] = np.exp(x[i])
  return ghost_config   


def frame_error(frame, ghostmap, blur_spectral=1, blur_spatial=1):
    fixed = fix_ghost_matrix(frame, ghostmap, blur_spectral = blur_spectral,
        blur_spatial=blur_spatial) 
    half = 640
    max_left = np.percentile(frame[:,:half],99)
    max_right = np.percentile(frame[:,half:],99)
    if max_left>max_right:
        return np.mean(pow(fixed[:,half:],2))# / np.mean(pow(frame[:,half:],2)) 
    else:
        return np.mean(pow(fixed[:,:half],2))# / np.mean(pow(frame[:,:half],2))


def err(x, frames, ghost_config, coarse):
    new_config = deserialize_ghost_config(x, ghost_config, coarse)
    ghostmap = build_ghost_matrix(new_config)
    blur_spatial = new_config['blur_spatial']
    blur_spectral = new_config['blur_spectral']
    jobs = [frame_error(frame, ghostmap, blur_spatial=blur_spatial,
         blur_spectral=blur_spectral) for frame in frames]
    errs = np.array(jobs)
    print(sum(errs))
    return sum(errs)
 
@ray.remote
def partial(x, i, frames, ghost_config, base_cost, coarse):
    x_perturb = x.copy()
    if coarse:
      eps = 0.01
    else:
      eps = 1e-7
    x_perturb[i] = x[i] + eps
    perturb_cost = err(x_perturb, frames, ghost_config, coarse)
    return (perturb_cost - base_cost)/eps


def jac(x, frames, ghost_config, coarse):
    base_cost = err(x,frames,ghost_config, coarse)
    jobs  = [partial.remote(x,i,frames,ghost_config,base_cost, coarse) for i in range(len(x))]
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
        frame = np.squeeze(I.load())
        if frame.shape[0] > frame.shape[1]:
            frame = frame.T
        frames.append(frame)
    frames = np.array(frames)
    with open(args.config,'r') as fin:
        ghost_config = json.load(fin)
 
    x0 = serialize_ghost_config(ghost_config, coarse=True)
    best = minimize(err, x0, args=(frames, ghost_config, True), jac=jac,method='SLSQP')
    best_config = deserialize_ghost_config(best.x, ghost_config, coarse=True)

    print(best.nit,'iterations')
    print('final error:',err(best.x, frames, ghost_config, coarse=True))
    print(best.message)

    x0 = serialize_ghost_config(best_config, coarse=False)
    best = minimize(err, x0, args=(frames, best_config, False), jac=jac,method='SLSQP')
    best_config = deserialize_ghost_config(best.x, best_config, coarse=False)

    print(best.nit,'iterations')
    print('final error:',err(best.x, frames, best_config, coarse=False))
    print(best.message)
    
    with open(args.output,'w') as fout:
        fout.write(json.dumps(best_config,indent=2))

if __name__ == '__main__':
    main()
