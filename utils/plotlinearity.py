# David R Thompson
import argparse, sys, os
from spectral.io import envi
import numpy as np
import pylab as plt

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


left, right, top, bottom = 25, 1264, 26, 313

def main():

    description = "Plot Linearity Curves"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    args = parser.parse_args()

    xs,ys = [],[]
    nfiles = len(args.input) 
    illums =[] 
    data = []

    for fi,infilepath in enumerate(args.input):

        infile = envi.open(find_header(infilepath))
        x = np.squeeze(infile.load())

        for i in range(0,x.shape[0]):
            plt.plot(x[i,:])
        plt.show()

if __name__ == '__main__':

    main()
