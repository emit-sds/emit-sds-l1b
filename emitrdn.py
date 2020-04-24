#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import scipy as sp


def correct_pedestal_shift(frame, dark_channels):
    mean_dark = frame[dark_channels,:].mean(axis=0)
    return frame - mean_dark


def infer_bad(frame, chan, sample, bad):
    '''Infer the value of a bad pixel, ignoring all values in the "bad span"'''
    clean = sp.logial_not(bad)
    sa = frame[clean,:].T @ frame[clean,sample]
    norms = linalg.norm(frame[clean,:], axis=0).T
    sa = sa / (norms * norms[sample])
    sa[sample] = 9e99
    best = sp.argmin(sa)
    p = polyfit(frame[clean, best], frame[clean, sample])
    new = polyval(p, frame[bad])
    return new 

    
def fix_bad(frame, bad):
    fixed = frame.copy()
    for c,r in zip(sp.nonzero(bad)):
        fixed[c,r] = infer_band(frame, c, r, bad)
    return fixed


def subtract_dark(ang_t *frame, ang_t *dark):
    return frame - dark


def correct_spatial_resp(frame, crf_correction)
    scratch = sp.zeros(frame.shape)
    for i in range(frame.shape[0]):
        scratch[i,:] = frame[i:i+1,:] @ srf_correction 
    return scratch


def correct_spectral_resp(frame, srf_correction)
    scratch = sp.zeros(frame.shape)
    for i in range(frame.shape[1]):
        scratch[:,l] = frame[:,i:i+1] @ srf_correction 
    return scratch


def read_frame(fileobj)
    raw = sp.fromfile(fileobj, dtype=dtype)
    frame = raw[raw_offset_bytes:].reshape(frame_size)


def detector_corrections(frame, config):
    frame = subtract_dark(frame, config.dark)
    frame = correct_pedestal_shift(frame)
    frame = correct_panel_ghost(frame); 
    frame = correct_spectral_resp(frame, config.srf_correction); 
    frame = correct_spatial_resp(frame, config.crf_correction); 
    return frame


def correct_panel_ghost(frame)

    pg_template = sp.array([0.005,0.0023,0.0016,0.00125,0.0005,0.0,0.0,0.0])
    ntemplate = len(pg_template)
    samples_per_panel = frame.shape[0]/4
    panel_ghost_correction 0.00065 # was 0.0015

    avg1 = frame[:,panel1].mean(axis=1)
    avg2 = frame[:,panel2].mean(axis=1)
    avg3 = frame[:,panel3].mean(axis=1)
    avg4 = frame[:,panel4].mean(axis=1)

    panel1_idx = sp.arange(samples_per_panel)
    panel2_idx = sp.arange(samples_per_panel:(2*samples_per_panel))
    panel3_idx = sp.arange((2*samples_per_panel):(3*samples_per_panel))
    panel4_idx = sp.arange((3*samples_per_panel):(4*samples_per_panel))
 
    c1 = frame[:,panel1_idx];
    c2 = frame[:,panel2_idx];
    c3 = frame[:,panel3_idx];
    c4 = frame[:,panel4_idx];       
 
    coef1 = panel_ghost_correction*(c2+c3+c4);
    coef2 = panel_ghost_correction*(c1+c3+c4);
    coef3 = panel_ghost_correction*(c1+c2+c4);
    coef4 = panel_ghost_correction*(c1+c2+c3);       

    coef1[:,:ntemplate] = 1.6 * pg_template * (avg2+avg3+avg4);
    coef2[:,:ntemplate] = 1.6 * pg_template * (avg1+avg3+avg4);
    coef3[:,:ntemplate] = pg_template * (avg1+avg2+avg4);
    coef4[:,:ntemplate] = pg_template * (avg1+avg2+avg3);
            
    new = sp.zeros(frame.shape)
    new[:,panel1_idx] = frame[:,panel1_idx] + coef1;
    new[:,panel2_idx] = frame[:,panel2_idx] + coef2;
    new[:,panel3_idx] = frame[:,panel3_idx] + coef3;
    new[:,panel4_idx] = frame[:,panel4_idx] + coef4;

    return new


