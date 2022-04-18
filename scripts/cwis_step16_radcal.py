# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from scipy.interpolate import interp1d

plot = True

q,wl,fwhm = np.loadtxt('../data/CWIS_Wavelengths_20220331.txt').T * 1000.0

# This helper function resamples a spectrum to our new wavelengths
def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)

# Load irradiance, translate uncertainty from percenmt to one sigma irradiance
irr = np.loadtxt('../data/ogse/cwis2/lamp_s1305_irradiance.txt', skiprows=1)
wl_irr = np.arange(250,2501)
irr = irr * 1e6 # translate W/(cm^2 nm) to microwatts /(cm^2 nm)
irradiance = resample(wl_irr, irr, method='cubic')

# Spectralon reflectance
wl_spec, spec_rfl = \
     np.loadtxt('../data/ogse/cwis2/panel_srt-99-120-0715_reflectance.txt',skiprows=1).T  
spectralon_rfl = resample(wl_spec, spec_rfl)

# Spectralon panel BRDF
brdf_factor = np.ones(len(wl)) * 1.015

# Radiance predicted at sensor
rdn = irradiance * spectralon_rfl / np.pi * brdf_factor 

# Wavelengths used for uncertainty data
wl_uncert = np.array([350.0, 450.0, 555.0, 654.6, 900.0, 1000.0, 2000.0, 2300.0, 2400.0])

# The NIST lamp uncertainty, courtesy Mark Helmlinger
nist_uncert_pct = resample(wl_uncert, np.array([1.27, 0.91, 0.77, 0.69, 0.57, 0.47, 0.50, 0.49, 1.11]))
nist_uncert_frac = nist_uncert_pct / 100.0 
nist_uncert = nist_uncert_frac * rdn # radiance units

# OGSE Geometry (conservative)
setup_uncert_pct = resample(wl_uncert, np.array([2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00 ]))
setup_uncert_frac = setup_uncert_pct / 100.0 
setup_uncert = setup_uncert_frac * rdn # radiance units

# Manufacturer transfer uncertainty (conservative)
xfer_uncert_pct = resample(wl_uncert, np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 ]))
xfer_uncert_frac = setup_uncert / 100.0
xfer_uncert = setup_uncert_frac * rdn # radiance units

# Lamp current uncertainty, from MH 
current_uncert_pct = resample(wl_uncert, np.array([0.71, 0.56, 0.46, 0.40, 0.30, 0.27, 0.14, 0.13, 0.12]))
current_uncert_frac = current_uncert_pct / 100.0
current_uncert = current_uncert_frac * rdn

# Spectral calibration-related uncertainy, from MH
spectral_uncert_pct = resample(wl_uncert, np.array([0.55, 0.22, 0.10, 0.06, 0.02, 0.01, 0.00, 0.00, 0.00]))
spectral_uncert_frac = spectral_uncert_pct / 100.0 
spectral_uncert = spectral_uncert_frac * rdn

# Spectralon reflectance uncertainty from manufaturer, reported as two standard deviations
spectralon_uncert_k2_frac = resample(wl_uncert, np.array([0.0055, 0.0055, 0.0055, 0.0049, 0.0049, 0.0049, 0.0088, 0.032, 0.032])) 
spectralon_uncert_frac = spectralon_uncert_k2_frac / 2.0 
spectralon_uncert = spectralon_uncert_frac * rdn

# BRDF uncertainty, conservative
brdf_uncert_pct = resample(wl_uncert, np.array([1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]))
brdf_uncert_frac = brdf_uncert_pct / 100.0
brdf_uncert = brdf_uncert_frac  * rdn

# Radiance uncertainty is Root Sum Square of independent sources
rdn_uncert = np.sqrt(nist_uncert**2 + setup_uncert**2 + \
        xfer_uncert**2 + spectral_uncert**2 + brdf_uncert**2 +\
        current_uncert**2 + spectralon_uncert**2)
rdn_uncert_frac = rdn_uncert / rdn

# Get the average DNs from each channel of the observation
basedir = '../data/'
input_file = basedir+'CWIS_FlatField_20220413'
I = envi.open(input_file+'.hdr')
DN = np.array([float(d) for d in I.metadata['average_dns']])
DN_std = np.array([float(d) for d in I.metadata['stdev_dns']])
channels = np.arange(len(wl),dtype=int)

# Calibration factors and their uncertainties
factors = rdn / DN
factors_uncert = rdn_uncert / DN

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
c = np.argmin(abs(wl-1440))
d = np.argmin(abs(wl-1450))
edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
interior = np.arange(c,b)
channels = np.arange(len(factors))
model = np.polyfit(channels[edges],factors[edges],2)
factors[interior] = np.polyval(model, channels[interior])

# Show the interpolated result
plt.plot(wl,factors)
plt.ylim([0,0.002])
plt.xlim([380,2500])
plt.show()

#factors_uncert = rdn_uncert / DN
SNR = DN/DN_std


if True:
    # These filenames are used for the automatic selection method
    timestamp = '20220413'
    np.savetxt('../data/CWIS_RadiometricCoeffs_'+timestamp+'.txt',
              np.c_[channels,factors,factors_uncert], fmt='%10.8f')
rdn_uncert = np.sqrt(nist_uncert**2 + setup_uncert**2 + \
        xfer_uncert**2 + spectral_uncert**2 + brdf_uncert**2 +\
        current_uncert**2 + spectralon_uncert**2)


    np.savetxt('../data/CWIS_RadiometricUncertainty_'+timestamp+'.csv',
              np.c_[wl,nist_uncert_frac,setup_uncert_frac,xfer_uncert_frac,
                      spectral_uncert_frac,brdf_uncert_frac,current_uncert_frac,
                      spectralon_uncert_frac,rdn_uncert_frac], fmt='%10.8f', delimiter=',',
              header='Wavelength, NIST calibration uncertainty, OGSE geometry uncertainty, Manufacturer transfer uncertainty, Spectral calibration uncertainty, BRDF uncertainty, Lamp current uncertainty, Spectralon reflectance uncertainty, Total uncertainty')
    np.savetxt('../data/CWIS_RadiometricReference_'+timestamp+'.txt',
              np.c_[wl,rdn], fmt='%10.8f')
    np.savetxt('../data/CWIS_RadiometricReferenceDN_'+timestamp+'.txt',
              np.c_[wl,DN], fmt='%10.8f')
    np.savetxt('../data/CWIS_RadiometricReferenceSNR_'+timestamp+'.txt',
              np.c_[wl,SNR], fmt='%10.8f')

if plot:
   plt.plot(wl,irradiance,color=[0.8,0.2,0.2])
   plt.plot(wl_irr,irr,'ko')
   plt.show()

if plot:
   plt.figure(figsize=(9,9))
   plt.plot(wl, nist_uncert / rdn * 100, color=[0.8, 0.2, 0.2])
   plt.plot(wl, spectralon_uncert / rdn * 100, color=[0.2, 0.8, 0.2])
   plt.plot(wl, setup_uncert / rdn * 100, color=[0.2, 0.2, 0.8])
   plt.plot(wl, current_uncert / rdn * 100, color=[0.2, 0.8, 0.8])
   plt.plot(wl, brdf_uncert / rdn * 100, color=[0.8, 0.2, 0.8])
   plt.plot(wl, spectral_uncert / rdn * 100, color=[0.8, 0.8, 0.2])
   plt.plot(wl, xfer_uncert / rdn * 100, color=[0.8, 0.8, 0.8])
   plt.plot(wl, rdn_uncert / rdn * 100, color='k')

   plt.legend(('Lamp calibration','Spectralon reflectance', 'OGSE geometry','Lamp current','Spectralon BRDF','Spectral calibration','Manufacturer transfer','Total uncertainty'))
   plt.grid(True)
   plt.box(False)
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Radiometric uncertainty (%)')
   plt.xlim([380,2500])
   plt.ylim([0,20])
   plt.show()


