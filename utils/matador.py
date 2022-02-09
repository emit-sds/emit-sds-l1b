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

colors = ['k','r','b']


def main():

    description = "Ratio two spectra"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    args = parser.parse_args()

    for i in range(0,len(args.input),2):

        infileA = envi.open(find_header(args.input[i]))
        infileB = envi.open(find_header(args.input[i+1]))
            
        A = infileA.load().mean(axis=0).T
        B = infileB.load().mean(axis=0).T
        print(B.shape,A.shape)
        
        acol = np.argmax(np.sum(A,axis=0))
        bcol = np.argmax(np.sum(B,axis=0))
        
        margin=5
        a = A[:,(acol-margin):(acol+margin+1)].mean(axis=1)
        b = B[:,(acol-margin):(acol+margin+1)].mean(axis=1)
        plt.plot(a/b,colors[i])

    plt.show()

if __name__ == '__main__':

    main()
