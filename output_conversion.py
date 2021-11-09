"""
This code contains support code for formatting L1B products for the LP DAAC.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
         Nimrod Carmon, nimrod.carmon@jpl.nasa.gov
"""

import argparse
from netCDF4 import Dataset
from emit_utils.daac_converter import add_variable, makeDims, makeGlobalAttr, add_loc, add_glt
from emit_utils.file_checks import netcdf_ext, envi_header
from spectral.io import envi
import os
import logging


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L1B PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('output_filename', type=str, help="Output netcdf filename")
    parser.add_argument('rdn_file', type=str, help="EMIT L1B radiance ENVI file")
    parser.add_argument('obs_file', type=str, help="EMIT L1B observation data ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
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
    loc_ds = envi.open(envi_header(args.loc_file))
    glt_ds = envi.open(envi_header(args.glt_file))

    # make the netCDF4 file
    logging.info(f'Creating netCDF4 file: {args.output_filename}')
    nc_ds = Dataset(args.output_filename, 'w', clobber=True, format='NETCDF4')

    # make global attributes
    logging.debug('Creating global attributes')
    makeGlobalAttr(nc_ds, args.rdn_file, args.glt_file)

    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis collection contains L1B at-sensor calibrated radiances \
        and geolocation data. The radiance calibration occurs in two basic stages: 1) transforming raw digital numbers \
        into radiance units using a calibrated radiometric response and correcting for electronic artifacts, and 2) \
        a correction for instrument optical effects which results in an absolute spectral wavelength calibration. \
        The geolocation data provides geotagging of the L1B radiance pixels (latitude, longitude, height) and \
        associated geometric information (path length, view angle, solar angle).\\n\\nThe calibrated at-sensor \
        radiance file contains radiance for each of {rdn_ds.shape} channels in units of microwatts per centimeter per \
        centimeter squared per steradian."
    nc_ds.sync()

    logging.debug('Creating dimensions')
    makeDims(nc_ds, args.rdn_file, args.glt_file)

    logging.debug('Creating and writing radiance metadata')
    add_variable(nc_ds, "sensor_band_parameters/radiance_wl", "f4", "Wavelength Centers", "nm",
                 [float(d) for d in rdn_ds.metadata['wavelength']], {"dimensions": ("number_of_bands",)})

    add_variable(nc_ds, "sensor_band_parameters/radiance_fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rdn_ds.metadata['fwhm']], {"dimensions": ("number_of_bands",)})

    #add_variable(nc_ds, "scan_line_attributes", "scan_start_time", "f8", "Scan start time (UTC)", "seconds",
    #             scan_start_times, {"dimensions": ("number_of_scans",)})
    #add_variable(nc_ds, "scan_line_attributes", "scan_end_time", "f8", "Scan start time (UTC)", "seconds",
    #             scan_end_times, {"dimensions": ("number_of_scans",)})

    #TODO: consider adding geotransform variable to each variable
    logging.debug('Creating and writing location data')
    add_loc(nc_ds, args.loc_file)
    logging.debug('Creating and writing glt data')
    add_glt(nc_ds, args.glt_file)

    logging.debug('Creating and writing obs data')
    nc_ds.createDimension('observation_bands', int(obs_ds.metadata['bands']))
    add_variable(nc_ds, "obs", "d", "Observation Data", None,
                 obs_ds.open_memmap(interleave='bip')[...].copy(), {"dimensions": ("ortho_y", "ortho_x", "observation_bands")})

    add_variable(nc_ds, "sensor_band_parameters/observation_bands", str, "Observation Band Names", None,
                 obs_ds.metadata['band names'], {"dimensions": ("observation_bands",)})

    add_variable(nc_ds, 'radiance', "f4", obs_ds.open_memmap(interleave='bip')[...].copy(), "Radiance Data",
                 "uW/cm/SR/nm", {"dimensions":("number_of_bands", "number_of_scans", "pixels_per_scan")})
    nc_ds.sync()
    nc_ds.close()
    logging.debug(f'Successfully created {args.output_filename}')

    return


if __name__ == '__main__':
    main()
