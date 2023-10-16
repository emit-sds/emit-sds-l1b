# David R Thompson
import numpy as np
import pylab as plt
import os, sys
from spectral.io import envi
from scipy.interpolate import interp1d

plot = True
frame_averaging = 15

q,wl,fwhm = np.loadtxt('../../data/cpm/CPM_Wavelengths_20231001.txt').T * 1000.0

def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)

# Window transmittance
wl_window, window_trans = np.loadtxt('../../data/ogse/cpm/window_infrasil301-302_transmittance.txt',skiprows=1).T
window_trans = resample(wl_window, window_trans)
window_uncert = np.ones(len(window_trans)) * 0.01


# -------------------------------------------------- NIST lamp and panel calibration

# Load irradiance, translate uncertainty from percenmt to one sigma irradiance
wl_irr, irr, irradiance_uncert = np.loadtxt('../../data/ogse/cpm/lamp_s1359_irradiance.txt', 
    skiprows=1).T 
irradiance = resample(wl_irr, irr, method='cubic')
irradiance_uncert = resample(wl_irr, irradiance_uncert)
irradiance_uncert = irradiance_uncert / 100.0 * irradiance

# Spectralon reflectance
wl_spec, spec_rfl, spec_uncert =\
     np.loadtxt('../../data/ogse/cpm/SRT-99-120-7015-2665.txt',skiprows=1).T  
spectralon_rfl = resample(wl_spec, spec_rfl)
spectralon_uncert = resample(wl_spec, spec_uncert)

# BRDF
brdf_factor = np.ones(len(wl)) * 0.995
brdf_uncert = np.ones(len(wl)) * 0.02

# integrated effect of r^-2 and cosine(i), uncerrtainty is 1% each  
fill_factor = np.ones(len(wl)) *  0.975 * 0.988  
fill_factor_uncert = np.ones(len(wl)) * 0.014
 
# Radiance 
rdn = irradiance * spectralon_rfl * window_trans / np.pi * brdf_factor * fill_factor
print('!',irradiance[10],spectralon_rfl[10],window_trans[10],brdf_factor[10])

distance_uncert = 0.0015875 # meters
distance = 0.5
distance_uncert_rdn =( 1-(0.5**2)/((0.5+distance_uncert)**2)) * rdn

# Derivatives of radiance
drdn_dirr    =              spectralon_rfl * window_trans / np.pi * brdf_factor * fill_factor 
drdn_dspec   = irradiance *                  window_trans / np.pi * brdf_factor * fill_factor 
drdn_dtrans  = irradiance * spectralon_rfl                / np.pi * brdf_factor * fill_factor 
drdn_dbrdf   = irradiance * spectralon_rfl * window_trans / np.pi               * fill_factor  
drdn_dfill   = irradiance * spectralon_rfl * window_trans / np.pi * brdf_factor                 

rdn_uncert = np.sqrt((drdn_dirr * irradiance_uncert)**2 + \
                     (drdn_dspec * spectralon_uncert)**2 + \
                     (drdn_dtrans * window_uncert)**2 + \
                     (drdn_dbrdf * brdf_uncert)**2 + \
                     (drdn_dfill * fill_factor_uncert)**2 + \
                     (distance_uncert_rdn**2))


basedir = '/beegfs/scratch/drt/20230928_CPM_Radcal/'

# subtract background
background = basedir+'20230904_193211_UTC_RADCAL/20230904_193416_UTC_RadCal_Fields-0-590_darksub_pedestal'
foreground = basedir+'20230904_193731_UTC_RADCAL/20230904_193936_UTC_RadCal_Fields-0-590_darksub_pedestal'
bgsub = basedir+'20230904_193731_UTC_RADCAL/20230904_193936_UTC_RadCal_Fields-0-590_darksub_pedestal_bgsub'
if not os.path.exists(bgsub):
  with open(background,'rb') as bg_fin:
    with open(foreground,'rb') as fg_fin:
      with open(bgsub,'wb') as fout:
        for i in range(14997):
          bg = np.fromfile(bg_fin, count=480*640, dtype=np.float32)
          fg = np.fromfile(fg_fin, count=480*640, dtype=np.float32)
          bgs = np.array(fg-bg,dtype=np.float32)
          bgs.tofile(fout)

# make header
bghdr = bgsub+'.hdr'
if not os.path.exists(bghdr):
  cmd = 'cp %s %s'%(foreground+'.hdr',bghdr)
  print(cmd)
  os.system(cmd)

# Run makeflat.py here
config = '../../config/cpm.json'          
flat = basedir+'20230904_193731_UTC_RADCAL/20230904_193936_UTC_RadCal_Fields-0-590_darksub_pedestal_bgsub_flat'
if True: #not os.path.exists(flat):
  cmd = 'python ../../utils/makeflat.py --config %s %s %s'%(config,bgsub,flat)
  print(cmd)
  os.system(cmd)

input_files = [flat]
output_files = ['CPM_RadiometricCoeffs_8ms_20231002.txt']

for input_file, outfile in zip(input_files, output_files):

    I = envi.open(input_file+'.hdr')
    DN = np.array([float(d) for d in I.metadata['average_dns']])
    DN_std = np.array([float(d) for d in I.metadata['stdev_dns']])

    channels = np.arange(len(wl),dtype=int)
    factors = rdn / DN

    plt.plot(wl,factors)

    # Interpolate over water vapor absorption at 1.88 microns
    a = np.argmin(abs(wl-1780))
    b = np.argmin(abs(wl-1800))
    c = np.argmin(abs(wl-1950))
    d = np.argmin(abs(wl-1970))
    edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
    interior = np.arange(c,b)
    channels = np.arange(len(factors))
    model = np.polyfit(channels[edges],factors[edges],2)
    factors[interior] = np.polyval(model, channels[interior])

    # Interpolate over water vapor absorption at 1.38 microns
    a = np.argmin(abs(wl-1335))
    b = np.argmin(abs(wl-1345))
    c = np.argmin(abs(wl-1420))
    d = np.argmin(abs(wl-1430))
    edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
    interior = np.arange(c,b)
    channels = np.arange(len(factors))
    model = np.polyfit(channels[edges],factors[edges],2)
    factors[interior] = np.polyval(model, channels[interior])

    # Show the interpolated result
    plt.semilogy(wl,factors)
    plt.ylim([0,0.02])
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Radiometric Calibration Coefficients (uW nm-1 sr-1 cm-2 DN-1)')
    plt.box(False)
    plt.grid(True)
    plt.xlim([380,2500])
    plt.show()

    factors_uncert = rdn_uncert / DN
    SNR = DN/DN_std/np.sqrt(frame_averaging)


    np.savetxt(outfile,np.c_[channels,factors,factors_uncert], fmt='%14.12f')

    if plot:
       plt.plot(wl,irradiance,color=[0.8,0.2,0.2])
       plt.plot(wl_irr,irr,'ko')
       plt.show()

    if plot:
       plt.figure(figsize=(9,9))
       plt.plot(wl, drdn_dspec * spectralon_uncert / rdn, color=[0.8, 0.2, 0.2])
       plt.plot(wl, drdn_dtrans * window_uncert / rdn, color=[0.2, 0.8, 0.2])
       plt.plot(wl, drdn_dirr * irradiance_uncert / rdn, color=[0.2, 0.2, 0.8])
       plt.plot(wl, drdn_dbrdf * brdf_uncert / rdn, color=[0.2, 0.8, 0.8])
       plt.plot(wl, drdn_dfill * fill_factor_uncert / rdn, color=[0.8, 0.8, 0.2])
       plt.plot(wl, distance_uncert_rdn / rdn, color=[0.8, 0.2, 0.8])
       plt.plot(wl, rdn_uncert/rdn, 'k--')
       plt.legend(('Spectralon reflectance','Window transmittance',
             'Lamp irradiance','Spectralon BRDF',
             'Pupil fill factor','OGSE geometry','Total uncertainty'))
       plt.grid(True)
       plt.box(False)
       plt.xlabel('Wavelength (nm)')
       plt.ylabel('Radiometric uncertainty, fractional')
       plt.xlim([380,2500])
       plt.ylim([0,0.1])
    plt.show()

cmd = 'cp ' 
cmd = cmd+flat
cmd = cmd + ' ' + '../../data/cpm/CPM_Flatfield_20231002'
print(cmd)
os.system(cmd)

cmd = 'cp ' 
cmd = cmd+flat+'.hdr'
cmd = cmd + ' ' + '../../data/cpm/CPM_Flatfield_20231002.hdr'
print(cmd)
os.system(cmd)

