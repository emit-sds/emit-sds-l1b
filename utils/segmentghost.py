from scipy.signal import find_peaks 
from scipy.stats import norm
from astropy import modeling
import sys
import numpy as np
from spectral.io import envi
from scipy.interpolate import interp1d
import json
from scipy.ndimage import gaussian_filter
import pylab as plt
from skimage.measure import LineModelND, ransac
np.set_printoptions(precision=5)

segment_points = [23,95,180,189,218,234,248,270,280]



with open(sys.argv[1],'r') as fin:
    config = json.load(fin)

new_config = {'blur':config['blur'],
               "center": config['center']}
for order in config['orders']:
    breaks = order['extent']
    slope = order['slope']
    offset = order['offset']
    for target_point in segment_points:
       source_point = (target_point - offset)/slope 
       if source_point > order['extent'][0] and source_point < order['extent'][1]:
           breaks.append(source_point)
    breaks.sort()
    for i in range(len(breaks-1)):
       new_config['orders'].append({'extent':(breaks[i],breaks[i]+1),
           'intensity':order['intensity'],
           'slope':slope,
           'offset':offset})

with open(sys.argv[2],'w' as fout:
    fout.write(json.dumps(new_config,indent=2))

