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
from numba import jit
from math import pow
from numba import jit


def fix_ghost(frame, fpa, ghostmap, blur_spatial, blur_spectral, center=649.5):

  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  if rows>cols:
      raise IndexError('Misformed frame')
  if blur_spatial>cols or blur_spectral>rows:
      raise IndexError('Illegal blur')

  for col in range(cols):
     tcol = int(center*2 - col)
     if tcol<0 or tcol>=fpa.native_columns:
         continue
     source = frame[:,col]
     target = source[np.newaxis,:] @ ghostmap
     ghost[:, tcol] = target 

  start = fpa.first_illuminated_row
  ghost[start:,:] = gaussian_filter(ghost[start:,:],[blur_spectral, blur_spatial])
  new = frame - ghost
  return new

