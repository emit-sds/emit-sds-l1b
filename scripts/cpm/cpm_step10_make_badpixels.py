# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
output_file = '../../data/cpm/CPM_BadElements_20231009'

rows, columns = 480, 640
bad_map = np.zeros((rows,columns,1), dtype=np.uint16)

envi.save_image(output_file+'.hdr',
    bad_map, interleave='bsq', ext='', force=True)
