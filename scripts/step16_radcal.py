# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from scipy.interpolate import interp1d

plot = True
frame_averaging = 15

q,wl,fwhm = np.loadtxt('../data/EMIT_Wavelengths_20211117.txt').T * 1000.0

def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)

# Load irradiance
# translate uncertainty from percenmt to one sigma irradiance
wl_irr, irr, irradiance_uncert = np.loadtxt('../data/ogse/tvac2/lamp_s1344_irradiance.txt', skiprows=2).T 
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

#plt.plot(wl,irradiance)
#plt.title('Irradiance')
#plt.show()
#plt.figure()
#plt.plot(wl,spectralon_rfl)
#plt.title('spectralon_rfl')
#plt.show()
#plt.figure()
#plt.plot(wl,mirror_rfl)
#plt.title('mirror_rfl')
#plt.show()
#plt.figure()
#plt.plot(wl,window_trans)
#plt.title('window_trans')
#plt.show()
#plt.figure()
#plt.plot(wl,brdf_factor)
#plt.title('brdf_factor')
#plt.show()

# Radiance 
rdn = irradiance * spectralon_rfl * mirror_rfl * window_trans / np.pi * brdf_factor

distance_uncert = 0.0015875 # meters
distance = 0.5
distance_uncert_rdn =( 1-(0.5**2)/((0.5+distance_uncert)**2)) * rdn

# Derivatives of radiance
drdn_dirr    =              spectralon_rfl * mirror_rfl * window_trans / np.pi * brdf_factor
drdn_dspec   = irradiance *                  mirror_rfl * window_trans / np.pi * brdf_factor
drdn_dtrans  = irradiance * spectralon_rfl * mirror_rfl                / np.pi * brdf_factor
drdn_dbrdf   = irradiance * spectralon_rfl * mirror_rfl * window_trans / np.pi 
drdn_dmirror = irradiance * spectralon_rfl *              window_trans / np.pi * brdf_factor

rdn_uncert = np.sqrt((drdn_dirr * irradiance_uncert)**2 + \
                     (drdn_dspec * spectralon_uncert)**2 + \
                     (drdn_dtrans * window_uncert)**2 + \
                     (drdn_dbrdf * brdf_uncert)**2 + \
                     (drdn_dmirror * mirror_uncert)**2 +\
                     (distance_uncert_rdn**2))

I = envi.open('../data/EMIT_FlatField_20211228.hdr')
DN = np.array([float(d) for d in I.metadata['average_dns']])
DN_std = np.array([float(d) for d in I.metadata['stdev_dns']])



channels = np.arange(len(wl),dtype=int)
factors = rdn / DN
factors_uncert = rdn_uncert / DN
SNR = DN/DN_std/np.sqrt(frame_averaging)
np.savetxt('../data/EMIT_RadiometricCoeffs_20211228.txt',
          np.c_[channels,factors,factors_uncert], fmt='%10.8f')
np.savetxt('../data/EMIT_RadiometricUncertainty_20211228.txt',
          np.c_[channels,factors_uncert/factors], fmt='%10.8f',
          header='Uncertainty, fractional')
np.savetxt('../data/EMIT_RadiometricReference_20211228.txt',
          np.c_[wl,rdn], fmt='%10.8f')
np.savetxt('../data/EMIT_RadiometricReferenceDN_20211228.txt',
          np.c_[wl,DN], fmt='%10.8f')
np.savetxt('../data/EMIT_RadiometricReferenceSNR_20211228.txt',
          np.c_[wl,SNR], fmt='%10.8f')



if plot:
   plt.plot(wl, drdn_dspec * spectralon_uncert / rdn, color=[0.8, 0.2, 0.2])
   plt.plot(wl, drdn_dtrans * window_uncert / rdn, color=[0.2, 0.8, 0.2])
   plt.plot(wl, drdn_dirr * irradiance_uncert / rdn, color=[0.2, 0.2, 0.8])
   plt.plot(wl, drdn_dbrdf * brdf_uncert / rdn, color=[0.2, 0.8, 0.8])
   plt.plot(wl, drdn_dmirror * mirror_uncert / rdn, color=[0.8, 0.8, 0.2])
   plt.plot(wl, distance_uncert_rdn / rdn, color=[0.8, 0.2, 0.8])
   plt.plot(wl, rdn_uncert/rdn, 'k')
   plt.legend(('Spectralon reflectance','Window transmittance',
         'Lamp irradiance','Spectralon BRDF','Mirror reflectance',
         'OGSE geometry','Total uncertainty'))
   plt.grid(True)
   plt.box(False)
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Radiometric uncertainty, fractional')
   plt.xlim([380,2500])
   plt.ylim([0,0.1])
   plt.show()
   #plt.plot(wl,DN/DN_std/np.sqrt(frame_averaging),'k')
   #plt.title('SNR')
   #plt.figure()
   #plt.plot(wl_irr,irr,'ko')
   #plt.plot(wl,irradiance,'r')
   #plt.figure()
   #plt.plot(wl_spec,spec_rfl,'ko')
   #plt.plot(wl,spectralon_rfl,'r')
   #plt.figure()
   #plt.plot(wl, factors)
   #plt.xlim([360,2540])
   #plt.figure()
   #plt.plot(wl[1:], np.diff(factors)/factors[1:])
   #plt.xlim([360,2540])
   #plt.figure()
   #plt.plot(wl, rdn_uncert/rdn)
   #plt.show()



