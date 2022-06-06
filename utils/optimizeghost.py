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
from fixghost import fix_ghost
from fixghostraster import build_ghost_matrix, build_ghost_blur
from fpa import FPA
import ray

rayargs={'num_cpus':40}
ray.init(**rayargs)


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def randomize_ghost_config(config, seed):
  new_config = deepcopy(config)
  if seed is None:
      return config
  rng = np.random.default_rng(int(seed))
  slopes = np.array([f['intensity_slope'] for f in config['orders']])
  offsets = np.array([f['intensity_offset'] for f in config['orders']])
  scalings = np.array([f['scaling'] for f in config['orders']])
  for i in range(len(config['orders'])):
     new_config['orders'][i]['scaling'] = \
       rng.random(None,float) * (scalings.max()-scalings.min()) + scalings.min()
     new_config['orders'][i]['intensity_offset'] = \
       rng.random(None,float) * (offsets.max()-offsets.min()) + offsets.min()
     new_config['orders'][i]['intensity_slope'] = \
       rng.random(None,float) * (slopes.max()-slopes.min()) + slopes.min()
  return new_config


def serialize_ghost_config(config, coarse):
  x, bounds = [],[]
  if coarse==1:
      for i in range(len(config['orders'])):
          x.append(config['orders'][i]['scaling'])
          bounds.append((0,9999))
  elif coarse==2:
      for i in range(len(config['orders'])):
          x.append(config['orders'][i]['intensity_slope'])
          bounds.append((-0.001,0.001))
          x.append(config['orders'][i]['intensity_offset'])
          bounds.append((-0.1,0.1))
  elif coarse==3:
      for i in range(len(config['orders'])):
          x.append(config['orders'][i]['extent'][0])
          bounds.append((x[-1]-2,x[-1]+2))
          x.append(config['orders'][i]['extent'][1])
          bounds.append((x[-1]-2,x[-1]+2))
  else:
      for zone in config['psf_zones']:
          for psf in zone['psfs']:
              x.append(psf['sigma'])
              bounds.append((0,100))
              x.append(psf['peak'])
              bounds.append((0,5))
  return x,bounds    


def deserialize_ghost_config(x, config, coarse):
  ghost_config = deepcopy(config) 
  if coarse==1: 
      for i in range(len(ghost_config['orders'])):
        ghost_config['orders'][i]['scaling'] = x[i]
  elif coarse==2: 
      for i in range(len(ghost_config['orders'])):
        ghost_config['orders'][i]['intensity_slope'] = x[i*2]
        ghost_config['orders'][i]['intensity_offset'] = x[i*2+1]
  elif coarse==3:
      for i in range(len(config['orders'])):
          ghost_config['orders'][i]['extent'][0] = x[i*2]
          ghost_config['orders'][i]['extent'][1] = x[i*2+1]
  else:
      ind = 0
      for zone in range(len(ghost_config['psf_zones'])):
          for psf in range(len(ghost_config['psf_zones'][zone]['psfs'])):
              ghost_config['psf_zones'][zone]['psfs'][psf]['sigma'] = x[ind]
              ghost_config['psf_zones'][zone]['psfs'][psf]['peak'] = x[ind+1]
              ind = ind+2
  return ghost_config   


def frame_error(frame, fpa, ghostmap, blur, center):
    try:
        fixed = fix_ghost(frame, fpa, ghostmap, blur=blur, center=center, plot=False) 
    except IndexError:
         # Something is out of bounds
         return 9e99
    half = int(round(fpa.native_columns/2))
    max_left = np.percentile(frame[:,:half],99)
    max_right = np.percentile(frame[:,half:],99)
    if max_left>max_right:
        return np.mean(pow(fixed[:,half:],2))# / np.mean(pow(frame[:,half:],2)) 
    else:
        return np.mean(pow(fixed[:,:half],2))# / np.mean(pow(frame[:,:half],2))


def err(x, fpa, frames, ghost_config, coarse):
    new_config = deserialize_ghost_config(x, ghost_config, coarse)
    ghostmap = build_ghost_matrix(new_config, fpa)
    blur = build_ghost_blur(new_config, fpa)
    center = new_config['center']
    jobs = [frame_error(frame, fpa, ghostmap, blur,
          center) for frame in frames]
    errs = np.array(jobs)
    print(sum(errs))
   #for i,err in enumerate(jobs):
   #    print('frame %i error %10.2f'%(i,err))
    return sum(errs)
 

@ray.remote
def partial(x, i, fpa, frames, ghost_config, base_cost, coarse):
    x_perturb = x.copy()
    if coarse:
        eps = 1e-7
    else:
        eps = 0.001
    x_perturb[i] = x[i] + eps
    perturb_cost = err(x_perturb, fpa, frames, ghost_config, coarse)
    return (perturb_cost - base_cost)/eps


def jac(x, fpa, frames, ghost_config, coarse):
    base_cost = err(x,fpa, frames,ghost_config, coarse)
    jobs  = [partial.remote(x,i,fpa,frames,ghost_config,base_cost, coarse) for i in range(len(x))]
    derivs = ray.get(jobs)
    return np.array(derivs)


def main():

    description = "Optimize ghost model"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('ghost_config')
    parser.add_argument('--config', default=None)
    parser.add_argument('--seed', default=None)
    parser.add_argument('--method', default='TNC')
    parser.add_argument('input',nargs='+')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    frames = []
    for infile in args.input:
        I = envi.open(find_header(infile))
        frame = np.squeeze(I.load())
        if frame.shape[0] > frame.shape[1]:
            frame = frame.T
        frames.append(frame)
    frames = np.array(frames)

    with open(args.ghost_config,'r') as fin:
        ghost_config = json.load(fin)

    ghost_config = randomize_ghost_config(ghost_config, args.seed)

    # We perform coordinate descent on different state vector subspaces
    for coarse in [3,1,2,0,1,2,0]:#[2,2]:#[1,0,2,1,0,2,1,0,2]:

        # nonlinear solution
        x0, bounds = serialize_ghost_config(ghost_config, coarse=coarse)
        best = minimize(err, x0, args=(fpa, frames, ghost_config, coarse), \
            jac=jac, bounds=bounds, method=args.method)
        best_config = deserialize_ghost_config(best.x, ghost_config, \
            coarse=coarse)
        
        # Print the result to screen
        print(best.nit,'iterations')
        print('final error:',err(best.x, fpa, frames, ghost_config, coarse=coarse))
        print(best.message)

        # Record final error
        xbest, bounds = serialize_ghost_config(best_config, coarse=2)
        best_config['final_error'] = err(xbest, fpa, frames, ghost_config, coarse=2)
        
        # Write provisional configuration to the output file
        with open(args.output,'w') as fout:
            fout.write(json.dumps(best_config,indent=2))

        # Initialize for the next round
        ghost_config = best_config

if __name__ == '__main__':
    main()
