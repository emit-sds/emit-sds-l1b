# David R Thompson
import os, sys, argparse
import numpy as np
import pylab as plt
from spectral.io import envi 


def main():

    description = "Calculate DN and Noise"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('DN')
    parser.add_argument('noise')
    args = parser.parse_args()

    infile = envi.open(args.input+'.hdr')
    data = infile.load() 

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

    maxrow = np.argmax(np.median(np.median(data, axis=2), axis=0))
    print('Max row',maxrow)
    maxchan = np.argmax(np.median(np.median(data, axis=1), axis=0))
    print('Max chan',maxchan)
    DN = np.mean(data[:,maxrow,:],axis=0)
    noise = np.std(data[:,maxrow,:],axis=0)
    np.savetxt(args.DN, DN, fmt='%10.8f')
    np.savetxt(args.noise, noise, fmt='%10.8f')

if __name__ == '__main__':

    main()
