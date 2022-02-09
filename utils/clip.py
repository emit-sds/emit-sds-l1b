# David R Thompson
import numpy as np
import os, os.path, sys
import json
import argparse
import logging
from spectral.io import envi

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def clip_frame(frame, start_row, end_row, start_col, end_col):
    clipped = frame[start_row:(end_row+1), :]
    clipped = clipped[:, start_col:(end_col+1)]
    return clipped


def main():

    description = "Fix pedestal shift for a data cube"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--start_row',default=None,type=int)
    parser.add_argument('--start_col',default=None,type=int)
    parser.add_argument('--end_row',default=None,type=int)
    parser.add_argument('--end_col',default=None,type=int)
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))

    if int(infile.metadata['data type']) == 2:
        dtype = np.int16
    elif int(infile.metadata['data type']) == 12:
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

    start_row = 0
    if args.start_row is not None:
        start_row = args.start_row
    start_col = 0
    if args.start_col is not None:
        start_col = args.start_col
    end_row = rows-1
    if args.end_row is not None:
        end_row = args.end_row
    end_col = columns-1
    if args.end_col is not None:
        end_col = args.end_col


    metadata = infile.metadata.copy()
    metadata['data type'] = 4
    metadata['bands'] = end_row - start_row + 1
    metadata['samples'] = end_col - start_col + 1
    envi.write_envi_header(args.output+'.hdr', metadata)

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        for line in range(lines):

            # Read a frame of data
            if line%10==0:
                logging.info('Line '+str(line))
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            clipped = clip_frame(frame, start_row, end_row, start_col, end_col)
            np.array(clipped, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
 

