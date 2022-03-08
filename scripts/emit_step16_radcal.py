# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from scipy.interpolate import interp1d

plot = True
frame_averaging = 15

q,wl,fwhm = np.loadtxt('../data/EMIT_Wavelengths_20220117.txt').T * 1000.0

def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)

# Load irradiance, translate uncertainty from percenmt to one sigma irradiance
wl_irr, irr, irradiance_uncert = np.loadtxt('../data/ogse/tvac4/lamp_s1352_irradiance.txt', skiprows=2).T 
irradiance = resample(wl_irr, irr, method='cubic')
irradiance_uncert = resample(wl_irr, irradiance_uncert)
irradiance_uncert = irradiance_uncert / 100.0 * irradiance

# Mirror transmittance
wl_mirror, mirror_rfl = np.loadtxt('../data/ogse/tvac2/mirror_coating_reflectance.txt').T
mirror_rfl = resample(wl_mirror, mirror_rfl)
mirror_uncert = np.ones(len(mirror_rfl)) * 0.01

# Spectralon reflectance
wl_spec, spec_rfl, spec_uncert =\
     np.loadtxt('../data/ogse/tvac2/panel_srt-99-120_reflectance.txt',skiprows=2).T  
spectralon_rfl = resample(wl_spec, spec_rfl)
spectralon_uncert = resample(wl_spec, spec_uncert)

# Window transmittance
wl_window, window_trans = np.loadtxt('../data/ogse/tvac2/window_infrasil301-302_transmittance.txt',skiprows=1).T
window_trans = resample(wl_window, window_trans)
window_uncert = np.ones(len(window_trans)) * 0.01

# BRDF
brdf_factor = np.ones(len(wl)) * 1.015
brdf_uncert = np.ones(len(wl)) * 0.01

# Transfer calibration
xfer_factor = 1.0
xfer_uncert = 0.0

# Radiance 
rdn = irradiance * spectralon_rfl * mirror_rfl * window_trans / np.pi * brdf_factor * xfer_factor
print('!',irradiance[10],spectralon_rfl[10],mirror_rfl[10],window_trans[10],brdf_factor[10])

distance_uncert = 0.0015875 # meters
distance = 0.5
distance_uncert_rdn =( 1-(0.5**2)/((0.5+distance_uncert)**2)) * rdn

# Derivatives of radiance
drdn_dirr    =              spectralon_rfl * mirror_rfl * window_trans / np.pi * brdf_factor * xfer_factor
drdn_dspec   = irradiance *                  mirror_rfl * window_trans / np.pi * brdf_factor * xfer_factor
drdn_dtrans  = irradiance * spectralon_rfl * mirror_rfl                / np.pi * brdf_factor * xfer_factor
drdn_dbrdf   = irradiance * spectralon_rfl * mirror_rfl * window_trans / np.pi               * xfer_factor 
drdn_dmirror = irradiance * spectralon_rfl *              window_trans / np.pi * brdf_factor * xfer_factor 
drdn_dxfer   = irradiance * spectralon_rfl * mirror_rfl * window_trans / np.pi * brdf_factor 

rdn_uncert = np.sqrt((drdn_dirr * irradiance_uncert)**2 + \
                     (drdn_dspec * spectralon_uncert)**2 + \
                     (drdn_dtrans * window_uncert)**2 + \
                     (drdn_dbrdf * brdf_uncert)**2 + \
                     (drdn_dmirror * mirror_uncert)**2 +\
                     (drdn_dxfer * xfer_uncert)**2 +\
                     (distance_uncert_rdn**2))

basedir = '/beegfs/scratch/drt/20220303_EMIT_TVAC4b/radcal/radcal/'
input_file = basedir+'emit20220305t002601_o00000_s000_l1a_raw_b0100_v01_strip_shift_darksub_pedestal_flat'
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
a = np.argmin(abs(wl-1340))
b = np.argmin(abs(wl-1350))
c = np.argmin(abs(wl-1410))
d = np.argmin(abs(wl-1420))
edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
interior = np.arange(c,b)
channels = np.arange(len(factors))
model = np.polyfit(channels[edges],factors[edges],2)
factors[interior] = np.polyval(model, channels[interior])

# Show the interpolated result
plt.plot(wl,factors)
plt.ylim([0,0.001])
plt.show()

factors_uncert = rdn_uncert / DN
SNR = DN/DN_std/np.sqrt(frame_averaging)


if True:
    # These filenames are used for the automatic selection method
    timestamp = '20220308'
    np.savetxt('../data/EMIT_RadiometricCoeffs_'+timestamp+'.txt',
              np.c_[channels,factors,factors_uncert], fmt='%10.8f')
    np.savetxt('../data/EMIT_RadiometricUncertainty_'+timestamp+'.txt',
              np.c_[channels,factors_uncert/factors], fmt='%10.8f',
              header='Uncertainty, fractional')
    np.savetxt('../data/EMIT_RadiometricReference_'+timestamp+'.txt',
              np.c_[wl,rdn], fmt='%10.8f')
    np.savetxt('../data/EMIT_RadiometricReferenceDN_'+timestamp+'.txt',
              np.c_[wl,DN], fmt='%10.8f')
    np.savetxt('../data/EMIT_RadiometricReferenceSNR_'+timestamp+'.txt',
              np.c_[wl,SNR], fmt='%10.8f')


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
   plt.plot(wl, drdn_dmirror * mirror_uncert / rdn, color=[0.8, 0.8, 0.2])
   plt.plot(wl, distance_uncert_rdn / rdn, color=[0.8, 0.2, 0.8])
   plt.plot(wl, drdn_dxfer * xfer_uncert / rdn,'k')
   plt.plot(wl, rdn_uncert/rdn, 'k--')
   plt.legend(('Spectralon reflectance','Window transmittance',
         'Lamp irradiance','Spectralon BRDF','Mirror reflectance',
         'OGSE geometry','Transfer uncertainty','Total uncertainty'))
   plt.grid(True)
   plt.box(False)
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Radiometric uncertainty, fractional')
   plt.xlim([380,2500])
   plt.ylim([0,0.1])
   plt.show()



