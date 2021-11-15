# David R Thompson 
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from scipy.interpolate import interp1d
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.optimize import minimize
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def sum_of_gaussians(x, mean1=0, amplitude1=1., sigma1=1.,
                        amplitude2=1., sigma2=1.,
                        amplitude3=1., sigma3=1.):
    return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) +
            amplitude2 * np.exp(-0.5 * ((x - mean1) / sigma2)**2) +
            amplitude3 * np.exp(-0.5 * ((x - mean1) / sigma3)**2))


def err(x,v,obs):
    mdl = sum_of_gaussians(v,x[0],np.abs(x[1]),np.abs(x[2]),
                           np.abs(x[3]),np.abs(x[4]),
                           np.abs(x[5]),np.abs(x[6]))
    obs[obs<1e-6] = 1e-6
    err = np.sum(pow(np.log(obs)-np.log(mdl),2))
    return err


def find_peak(x):
    
    fitter = modeling.fitting.LevMarLSQFitter()

    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)   # depending on the data you need to give some initial values
    fitted_model = fitter(model, np.arange(len(x)), x)
    return fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]


def find_scatter(obs,plot=False):
    v = np.arange(len(obs))
    ctr = np.argmax(obs)
    x0 = np.array([ctr, 1, 0.7, 0.006, 2.5, 0.0007, 5])
    mdl = sum_of_gaussians(v,x0[0],abs(x0[1]),abs(x0[2]),abs(x0[3]),abs(x0[4]),abs(x0[5]),abs(x0[6]))
    xbest = minimize(lambda q: err(q, v, obs), x0, method='CG')
    x = xbest.x
    if plot:
        plt.semilogy(v,mdl,'b')
        plt.semilogy(v,obs,'ko')
        mdl = sum_of_gaussians(v,x[0],abs(x[1]),abs(x[2]),abs(x[3]),abs(x[4]),abs(x[5]),abs(x[6]))
        plt.semilogy(v,mdl,'r')
        plt.show()
    return x


def main():

    description = "Calculate Scatter function"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--spatial',action='store_true')
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--target_column',type=int,default=40)
    parser.add_argument('input',nargs='+')
    args = parser.parse_args()

    for infile in args.input:
        infile = envi.open(find_header(infile))
        
        if int(infile.metadata['data type']) == 2:
            dtype = np.uint16
        elif int(infile.metadata['data type']) == 4:
            dtype = np.float32
        else:
            raise ValueError('Unsupported data type')
        if infile.metadata['interleave'] != 'bil':
            raise ValueError('Unsupported interleave')
        
        rows = int(infile.metadata['bands'])
        columns = int(infile.metadata['samples'])
        lines = int(infile.metadata['lines'])
        nframe = rows * columns
        margin=20
        
        X = infile.load()
        c = np.argmax(np.sum(np.sum(X,axis=1),axis=0)) 
        
        # spatial scatter
        col = args.target_column
        sequence = X[margin:,(col-10):(col+11),:]
        sequence = np.mean(sequence, axis=2)
        sequence = np.mean(sequence, axis=0)
        sequence = np.squeeze(sequence)
        sequence = sequence / max(sequence)                                         
        best = find_scatter(sequence,args.plot)
        print('%i %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f'%(c, best[0],best[1],best[2],best[3],best[4],best[5],best[6])) 

if __name__ == '__main__':

    main()
