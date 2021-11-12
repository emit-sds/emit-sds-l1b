# David R Thompson
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.linalg import solve, inv
from scipy.stats import norm
import json
from numpy import convolve
from scipy.optimize import minimize
import optimparallel

frames = envi.open('frames.hdr').load()
frames = frames.transpose((0,2,1))
left, right, short, long = 25, 1265, 21, 314

with open('emit_ghost.json','r') as fin:
    ghost_config = json.load(fin)
    
def ghost_config_to_vector(ghost_config):
    x = [ghost_config['center']]
    for j,order in enumerate(ghost_config['orders']):
        x.append(order['extent'][0])
        x.append(order['extent'][1])
        x.append(order['slope'])
        x.append(order['offset'])
        x.append(order['blur_spatial'])
        x.append(order['blur_spectral'])
        x.append(order['intensity'])
    return np.array(x,dtype=float)

def vector_to_ghost_config(x):
    ghost_config = {}
    ghost_config['center'] = x[0]
    ghost_config['orders'] = []
    if (len(x)-1)%7 != 0:
        raise ValueError('Incorrect vector size')
    n_orders = int((len(x)-1)/7.0)
    for i in range(n_orders):
        order = {'extent': (x[i*7+1],x[i*7+2]),
                 'slope': x[i*7+3],
                 'offset': x[i*7+4],
                 'blur_spatial': x[i*7+5],
                 'blur_spectral': x[i*7+6],
                 'intensity': x[i*7+7]}
        ghost_config['orders'].append(order)
    return ghost_config

def get_ghost(frame, ghost_config):
    ghosted = np.zeros((480,1280))
    spectral_map = np.zeros((480,480))
    spatial_blurs = np.zeros(480)
    center = ghost_config['center']
    for i in range(480):
        for j,order in enumerate(ghost_config['orders']):
            extent = order['extent']
            slope, offset = order['slope'], order['offset']
            blur_spatial = order['blur_spatial']
            blur_spectral = order['blur_spectral']
            intensity = order['intensity']
            if i>=extent[0] and i<extent[1]:
                target = slope * i + offset
                spectral_map[i,:] = norm.pdf(range(480),target,blur_spectral)
                spectral_map[i,:] = spectral_map[i,:] / spectral_map[i,:].sum() * intensity 
                spatial_blurs[i] = blur_spatial
    wavelength_resorted = spectral_map.T @ frame 
    for col in range(1280):
        tcol = int(round(center+(center-col)))
        if tcol<20 or tcol>=1270:
            continue
        ghosted[:,tcol] = wavelength_resorted[:,col]
    for row in range(480):
        kernel = norm.pdf(np.arange(-4,5),0,spatial_blurs[row])
        kernel = kernel/kernel.sum()
        ghosted[:,tcol] = convolve(ghosted[:,tcol],kernel,mode='same')
    return ghosted
                                                            
def err(x, frames):
    half = 640
    config = vector_to_ghost_config(x)
    err = 0
    for i,integration in enumerate(frames):
        if i%25==0:
            print('error frame',i)
        frame = np.squeeze(integration)
        max_left = frame[:,:half].max()
        max_right = frame[:,half:].max()
        ghost = get_ghost(frame,config)
        if max_left>max_right:
            err = err + np.mean(pow(ghost[:,half:] - frame[:,half:],2))
        else:
            err = err + np.mean(pow(ghost[:,:half] - frame[:,:half],2))
    print(err)
    return err
  

x0 = ghost_config_to_vector(ghost_config)
best = minimize(lambda v: err(v,frames), x0)
best_config = vector_to_ghost_config(best.x)
  
with open('emit_ghost_optimized.json','w') as fout:
    json.dump(best_config,fout)           
                               
if False:    
    new = frame-ghost
    plt.figure(figsize=(15,15))
    plt.imshow(frame,vmin=-0.001,vmax=0.001)
    plt.figure(figsize=(15,15))
    plt.imshow(new,vmin=-0.001,vmax=0.001) 
    plt.figure(figsize=(15,15))
    plt.imshow(frame-new,vmin=-0.001,vmax=0.001)
    break
