# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from skimage.filters import threshold_otsu, difference_of_gaussians
from PIL import Image
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

# bad row, column ranges (inclusive, zero indexed, 480x1280 format]
manual_bads=[(60,164),
             (180,140), 
             (181,(140,141)), 
             (204,141),
             (351,73),
             (433,253), 
             (462,240),
             (488,222),
             (535,101),
             (740,266),
             (736,277),
             (811,202),
             (817,264),
             (911,51),
             (918,192),
             (988,218),
             (1054,37),
             (1172,129),
             (1205,271),
             (1232,81),
             (1258,131)]


input_file = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_80pc_darksub_pedestal'
mask_image = '/beegfs/scratch/drt/20220112_CWIS2/20220107_flatfield/light_80pc_mask.png'
output_file = '../data/CWIS_BadElements_20220413'
cue_channel=148


infile = envi.open(find_header(input_file))

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
margin=2

flat  = np.zeros((rows,columns))
count = np.zeros((rows,columns))
sumsq = np.zeros((rows,columns))
ref = np.zeros((lines,columns))
allctrs,alllines = [],[]

mask = np.asarray(Image.open(mask_image))
print(mask.shape) 
if len(mask.shape)>2:
    mask = np.array(mask[:,:,0]>0, dtype=int)
if mask.shape[0] != lines or mask.shape[1] != columns:
    raise IndexError('mask does not match image')


with open(input_file,'rb') as fin:

    for line in range(lines):

        # Read a frame of data
        frame = np.fromfile(fin, count=nframe, dtype=dtype)
        frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
        reference = frame[cue_channel, :]
        use = np.where(mask[line,:]>0)[0]

        if len(use)>0:
            print(line,len(use),np.median(use))#,thresh)
            flat[:,use] = flat[:,use] + frame[:,use] 
            count[:,use] = count[:,use] + 1
            sumsq[:,use] = sumsq[:,use] + pow(frame[:,use],2)

    mean_sumsq = sumsq / count
    flat = flat / count

    rowmean = flat[:,30:1250].mean(axis=1)
    rowstdev = flat[:,30:1250].std(axis=1)
    stdev = np.sqrt(mean_sumsq - pow(flat,2))
    stdev[np.logical_not(np.isfinite(stdev))] = 0

    stdev = difference_of_gaussians(stdev,0,2)
    plt.imshow(stdev,vmin=-100,vmax=100)
    plt.colorbar()
    plt.show()
    bad = stdev>90

bad[:,:23] = 0
bad[:,1265:] = 0

for bad_cols, bad_rows in manual_bads:
    if type(bad_rows)==int:
        rows_range = [bad_rows]
    else:
        rows_range = range(bad_rows[0],bad_rows[1]+1)
    if type(bad_cols)==int:
        cols_range = [bad_cols]
    else:
        cols_range = range(bad_cols[0],bad_cols[1]+1)
    for col in cols_range:
        for row in rows_range:
            if row>0 and row<bad.shape[0] and col>0 and col<bad.shape[1]:
                bad[row,col] = 1
            else:
                print('row ',row,' column ',col,' out of bounds')

if False:
    plt.hist(stdev.flatten(),500)
    plt.figure()
    plt.imshow(bad)
    plt.show()

bads = 0
bad_map = bad.copy()
bad_map = np.array(bad_map,dtype=np.int16)
for column in range(bad_map.shape[1]):
    state_machine = 0
    for row in range(bad_map.shape[0]):
        if bad[row,column]:
            state_machine = state_machine + 1
            bad_map[row,column] = -state_machine
            print(row,column,state_machine)
            bads = bads + 1
        else:
            state_machine = 0

print('total bads:',bads)
bad_map = bad_map.reshape((rows,columns,1))

envi.save_image(output_file+'.hdr',
    bad_map, interleave='bsq', ext='', force=True)
