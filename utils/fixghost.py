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
from emit_fpa import native_rows, frame_embed, frame_extract
from emit_fpa import first_illuminated_row
from numba import jit


@jit
def fix_ghost(frame, ghostmap, center=649.5, blur_spatial=50, blur_spectral=1, fudge = 4.25):

  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  if rows>cols:
      raise IndexError('Misformed frame')

  for row in range(rows):
    for ghost_row, intensity in ghostmap[row]:
       for col in range(cols):
          tcol = int(center*2 - col)
          if tcol>0 and tcol<1280:
              ghost[ghost_row, tcol] = \
                 ghost[ghost_row, tcol] + frame[row,col] * intensity * fudge

  start = first_illuminated_row
  ghost[start:,:] = gaussian_filter(ghost[start:,:],[blur_spectral, blur_spatial])
  new = frame - ghost
  return new

#@jit
def fix_ghost_matrix(frame, ghostmap, center=649.5, blur_spatial=50, blur_spectral=1, fudge = 4.25):

  ghost = np.zeros(frame.shape)
  rows, cols = frame.shape
  if rows>cols:
      raise IndexError('Misformed frame')

  for col in range(cols):
     tcol = int(center*2 - col)
     if tcol<0 or tcol>1279:
         continue
     source = frame[:,col]
     target = source[np.newaxis,:] @ ghostmap
     ghost[:, tcol] = target * fudge

  start = first_illuminated_row
  ghost[start:,:] = gaussian_filter(ghost[start:,:],[blur_spectral, blur_spatial])
  new = frame - ghost
  return new

