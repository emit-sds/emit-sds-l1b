#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: Philip G Brodrick, philip.brodrick@jpl.nasa.gov

import os, sys
import numpy as np
import logging
from numba import jit


#frame = fix_electronic_ghost(frame, fpa.eghost_samples_per_panel, np.array(fpa.eghost_template),
#                                     fpa.eghost_panel_correction, fpa.eghost_panel_multipliers)
    
def fix_electronic_ghost(frame, samples_per_panel, pg_template, panel_correction, panel_multipliers):

    if frame.shape[1] % samples_per_panel != 0:
        logging.error('Incorrect samples per panel')

    num_panels = int(frame.shape[1]/samples_per_panel)

    panel_means = [np.mean(frame[:,s*samples_per_panel:(s+1)*samples_per_panel],axis=1) for s in range(num_panels)]
    
    coefficients_means = [(np.sum(panel_means,axis=0) - panel_means[s]) for s in range(num_panels)] 

    across_panel_sum = np.sum(frame.reshape((frame.shape[0],-1,samples_per_panel)),axis=1) 
    
    coefficients_frame = np.zeros(frame.shape)
    for s in range(num_panels):
        coefficients_frame[:,s*samples_per_panel:(s+1)*samples_per_panel] = panel_correction*(across_panel_sum-frame[:,s*samples_per_panel:(s+1)*samples_per_panel])
        #Overwrite the first 8
        coefficients_frame[:,s*samples_per_panel:s*samples_per_panel+len(pg_template)] = panel_multipliers[s]*pg_template[np.newaxis,:]*coefficients_means[s][:,np.newaxis]
                                                                                              
    frame += coefficients_frame
    
    return frame




@jit(nopython=False)
def fix_electronic_ghost_dep(frame, samples_per_panel, pg_template, panel_correction, panel_multipliers):

    if frame.shape[0] % samples_per_panel != 0:
        logging.error('Incorrect samples per panel')
 
    num_panels = int(frame.shape[0]/samples_per_panel)
 
    panel_means = [np.mean(frame[:,s*samples_per_panel:(s+1)*samples_per_panel],axis=1) for s in range(num_panels)]
    coefficient_means = [(np.sum(panel_means,axis=0) - panel_means[s]) for s in range(num_panels)]

    fixed = frame.copy()

    for s in range(samples_per_panel):
        if s < len(pg_template):
            for n in range(num_panels):
                fixed[:, s + n*samples_per_panel] += panel_multipliers[n]*pg_template[s]*coefficient_means[n]

        else:
            c = [fixed[:,s + n*samples_per_panel] for n in range(num_panels)]
            coeff = [panel_correction * (np.sum(c,axis=0) - c[n]) for n in range(num_panels)]

            for n in range(num_panels):
                fixed[:, s + n*samples_per_panel] += coeff[n]

    return fixed



