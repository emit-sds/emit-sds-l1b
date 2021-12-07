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
import optimparallel

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


def serialize_ghost_config(config, expressive=False):

  x = [config['blur']]
  for i in range(len(config['orders'])):
      x.append(config['orders'][i]['extent'][0])
      x.append(config['orders'][i]['extent'][1])
      x.append(config['orders'][i]['intensity'])
      if expressive:
          if 'intensity_shift' in config['orders'][i]:
              for r in config['orders'][i]['intensity_shift']:
                  x.append(r)
          else:
              for j in range(int(round(config['orders'][i]['extent'][0])),
                             int(round(config['orders'][i]['extent'][1])+1)):
              
                  x.append(1.0)
  return x    


def deserialize_ghost_config(x, config, expressive=False):
  ghost_config = deepcopy(config) 
  if (not expressive) and (len(x)-1)/3 != len(config['orders']):
    raise IndexError('bad state vector size')
  ghost_config['blur'] = x[0]
  ind = 1
  for i in range(len(config['orders'])):
    ghost_config['orders'][i]['extent'][0] = int(round(x[ind]))
    ind = ind+1
    ghost_config['orders'][i]['extent'][1] = int(round(x[ind]))
    ind = ind+1
    ghost_config['orders'][i]['intensity'] = x[ind]
    ind = ind+1
    if expressive:
        ghost_config['orders'][i]['intensity_shift'] = []
        for j in range(int(round(ghost_config['orders'][i]['extent'][0])),
                       int(round(ghost_config['orders'][i]['extent'][1]+1))):
           ghost_config['orders'][i]['intensity_shift'].append(x[ind])
           ind = ind+1
           
  return ghost_config   


#@jit   
def fix_ghost(frame, config, expressive=False):

  center = config['center']
  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  blur = config['blur']

  for row in range(rows):
    for order in config['orders']:
       if row>=order['extent'][0] and row<=order['extent'][1]: 
          ghost_position = int(order['slope']*row + order['offset'])
          if ghost_position > 0 and ghost_position < 480:
              intensity = order['intensity']
              if expressive:
                 intensity_shift = order['intensity_shift'][int(round(row-order['extent'][0]))]
                 intensity = intensity * intensity_shift
              for col in range(cols):
                 tcol = int(center*2 - col)
                 if tcol>0 and tcol<1280:
                     ghost[ghost_position, tcol] = frame[row,col] * intensity

  ghost = gaussian_filter(ghost,[blur,blur])
  new = frame - ghost
  return new

@ray.remote
def frame_error(frame, new_config, expressive=False):
    fixed = fix_ghost(frame, new_config, expressive)
    half = 640
    max_left = frame[:,:half].max()
    max_right = frame[:,half:].max()
    if max_left>max_right:
        return np.mean(pow(fixed[:,half:],2))
    else:
        return np.mean(pow(fixed[:,:half],2))

#@ray.remote
def cost_gradient(x,i,frames,config,expressive=False):
    x_perturb = x.copy()
    x_perturb[i] = x_perturb[i]+eps
    print(len(x),len(x_purturb))
    cost = err(x, frames, config, expressive)
    cost_perturb = err(x_perturb, frames, config, expressive)
    return (cost_perturb-cost)/eps


def jac(x, frames, config, expressive=False):
    config = deserialize_ghost_config(x, config, expressive)
    jobs = [cost_gradient(x,i,frames,config,expressive=False) for i in range(len(x))] 
    jac = ray.get(jobs)
    return jac

def err(x, frames, config, expressive=False):
    new_config = deserialize_ghost_config(x, config, expressive)
    jobs = [frame_error.remote(frame, new_config, expressive) for frame in frames]
    errs = ray.get(jobs)
    if False:
         jobs = [frame_error(frame, new_config, expressive) for frame in frames]
         errs = np.array(jobs)
    print(sum(errs))
    return sum(errs)
  

def main():

    description = "Optimize ghost model"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config')
    parser.add_argument('input',nargs='+')
    parser.add_argument('--expressive',action='store_true')
    parser.add_argument('output')
    args = parser.parse_args()

    frames = []
    for infile in args.input:
        print(infile)
        I = envi.open(infile)
        frames.append(np.squeeze(I.load()))
    frames = np.array(frames)
    print(frames.shape)
    with open(args.config,'r') as fin:
        ghost_config = json.load(fin)

    x0 = serialize_ghost_config(ghost_config,expressive=args.expressive)
    opts = {'max_iters':10}
    best = minimize(err, x0, args=(frames, ghost_config, args.expressive), options=opts)
   #opts = {'max_workers':40}
   #opts = {'max_workers':20}
   #best = optimparallel.minimize_parallel(lambda v: err(v, frames, ghost_config, args.expressive), x0, parallel=opts)
    best_config = deserialize_ghost_config(best.x, ghost_config, args.expressive)
    
    with open(args.output,'w') as fout:
        fout.write(json.dumps(best_config,indent=2))

if __name__ == '__main__':
    main()
