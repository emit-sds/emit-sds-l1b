#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import json
import logging
import argparse
from scipy.ndimage import gaussian_filter
from math import pow
from pylab import plt


def fix_ghost(frame, fpa, ghostmap, blur, center, plot=False):

  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  if rows>cols:
      raise IndexError('Misformed frame')

  for col in range(cols):
     tcol = int(center*2 - col)
     if tcol<0 or tcol>=fpa.native_columns:
         continue
     source = frame[:,col]
     target = source[np.newaxis,:] @ ghostmap
     ghost[:, tcol] = target 

  for extent, (blur_spatial, blur_spectral) in blur.items():
      lo = max(extent[0],fpa.first_valid_row)
      hi = min(extent[1],fpa.last_valid_row)
      indices = np.arange(lo, hi+1)
      ghost = blur_spectral @ ghost
      ghost[indices,:] = (blur_spatial @ ghost[indices,:].T).T
  
  new = frame - ghost

  if plot:
      plt.imshow(frame,vmin=-10,vmax=50)
      plt.figure()
      plt.imshow(ghost,vmin=-10,vmax=50)
      plt.figure()
      plt.imshow(new,vmin=-10,vmax=50)
      plt.show()

  return new

