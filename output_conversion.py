"""
This code contains support code for formatting L1B products for the LP DAAC.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
         Nimrod Carmon, nimrod.carmon@jpl.nasa.gov
"""

import argparse
from netCDF4 import Dataset
from emit_utils import daac_converter 
from emit_utils.file_checks import netcdf_ext, envi_header
from spectral.io import envi
import os
import logging
import numpy as np


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L1B PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('rdn_output_filename', type=str, help="Radiance output netcdf filename")
    parser.add_argument('obs_output_filename', type=str, help="Observation output netcdf filename")
    parser.add_argument('rdn_file', type=str, help="EMIT L1B radiance ENVI file")
    parser.add_argument('obs_file', type=str, help="EMIT L1B observation data ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
    parser.add_argument('version', type=str, help="3 digit (with leading V) version number")
    parser.add_argument('software_delivery_version', type=str, help="The extended build number at delivery time")
    parser.add_argument('--flat_field_update', type=str, help="Optional flat field update")
    parser.add_argument('--ummg_file', type=str, help="Output UMMG filename")
    parser.add_argument('--log_file', type=str, default=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.log_level, filename=args.log_file)

    rdn_ds = envi.open(envi_header(args.rdn_file))
    obs_ds = envi.open(envi_header(args.obs_file))

    # Setup complete - transition to making the radiance output file

    logging.info(f'Creating Radiance netCDF4 file: {args.rdn_output_filename}')
    nc_ds = Dataset(args.rdn_output_filename, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    if args.flat_field_update:
        # Pass in runconfig in order to get ffupdate_paths
        runconfig_file = args.rdn_file.replace(".img", "_runconfig.json")
        daac_converter.makeGlobalAttr(nc_ds, args.rdn_file, args.software_delivery_version,
                                      glt_envi_file=args.glt_file, rdn_runconfig_file=runconfig_file)
    else:
        daac_converter.makeGlobalAttr(nc_ds, args.rdn_file, args.software_delivery_version,
                                      glt_envi_file=args.glt_file)

    nc_ds.title = "EMIT L1B At-Sensor Calibrated Radiance Data 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis file contains L1B at-sensor calibrated radiances. \
The radiance calibration occurs in two basic stages: 1) transforming raw digital numbers \
into radiance units using a calibrated radiometric response and correcting for electronic artifacts, and 2) \
a correction for instrument optical effects which results in an absolute spectral wavelength calibration. \
The radiance file contains radiance for each of {rdn_ds.shape[-1]} channels in units of microwatts per centimeter per \
centimeter squared per steradian. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    daac_converter.makeDims(nc_ds, args.rdn_file, args.glt_file)

    logging.debug('Creating and writing radiance metadata')
    daac_converter.add_variable(nc_ds, "sensor_band_parameters/wavelengths", "f4", "Wavelength Centers", "nm",
                 [float(d) for d in rdn_ds.metadata['wavelength']], {"dimensions": ("bands",)})
    daac_converter.add_variable(nc_ds, "sensor_band_parameters/fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rdn_ds.metadata['fwhm']], {"dimensions": ("bands",)})
    daac_converter.add_variable(nc_ds, 'radiance', "f4", "Radiance Data", "uW/cm^2/SR/nm", rdn_ds.open_memmap(interleave='bip')[...].copy(),
                 {"dimensions":("downtrack", "crosstrack", "bands")})
    if args.flat_field_update:
        ff_ds = envi.open(envi_header(args.flat_field_update))
        daac_converter.add_variable(nc_ds, 'flat_field_update', "f4", "Flat Field Update", None,
                                    ff_ds.open_memmap(interleave='bip')[...].copy().squeeze().T,
                                    {"dimensions": ("crosstrack", "bands")})

    logging.debug('Creating and writing location data')
    daac_converter.add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    daac_converter.add_glt(nc_ds, args.glt_file)

    nc_ds.sync()
    nc_ds.close()
    del nc_ds
    logging.debug(f'Successfully created {args.rdn_output_filename}')

    # Transition to creating the Observation output

    logging.info(f'Creating Observation netCDF4 file: {args.obs_output_filename}')
    nc_ds = Dataset(args.obs_output_filename, 'w', clobber=True, format='NETCDF4')

    logging.debug('Creating global attributes')
    daac_converter.makeGlobalAttr(nc_ds, args.obs_file, args.software_delivery_version, glt_envi_file=args.glt_file)

    nc_ds.title = "EMIT L1B Observation Data 60 m " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis file contains L1B geometric information (path length, view and solar angles, timing) associated with \
each pixel in an acquisition. \
Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included."
    nc_ds.sync()

    daac_converter.makeDims(nc_ds, args.obs_file, args.glt_file)


    logging.debug('Creating and writing obs data')
    daac_converter.add_variable(nc_ds, "obs", "f4", "Observation Data", None,
                 obs_ds.open_memmap(interleave='bip')[...].copy(), {"dimensions": ("downtrack", "crosstrack", "bands")})
    daac_converter.add_variable(nc_ds, "sensor_band_parameters/observation_bands", str, "Observation Band Names", None,
                 obs_ds.metadata['band names'], {"dimensions": ("bands",)})

    logging.debug('Creating and writing location data')
    daac_converter.add_loc(nc_ds, args.loc_file)

    logging.debug('Creating and writing glt data')
    daac_converter.add_glt(nc_ds, args.glt_file)

    nc_ds.sync()
    nc_ds.close()
    del nc_ds
    logging.debug(f'Successfully created {args.obs_output_filename}')


    return


if __name__ == '__main__':
    main()
