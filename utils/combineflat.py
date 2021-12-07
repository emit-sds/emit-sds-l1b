# David R Thompson
import argparse
from spectral.io import envi
import numpy as np
import sys, os


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')




def main():

    description = "Average and stack multiple images"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input',nargs='+')
    parser.add_argument('output')
    args = parser.parse_args()

    frames = []
    for fi,infilepath in enumerate(args.input):
        print(fi,'/',len(args.input))

        x = envi.open(infilepath+'.hdr').load()      
        x = np.squeeze(x)
        x = np.mean(axis=0)
        frames.append(x)

    frames = np.array(frames,dtype=np.float32)
    envi.save_image(args.output+'.hdr',frames,ext='',force=True)

if __name__ == '__main__':

    main()
