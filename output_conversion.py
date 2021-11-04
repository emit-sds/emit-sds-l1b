"""
This code contains support code for formatting L1B products for the LP DAAC.

Authors: Philip G. Brodrick, philip.brodrick@jpl.nasa.gov
         Nimrod Carmon, nimrod.carmon@jpl.nasa.gov
"""

import argparse
from netCDF4 import Dataset
from emit_utils.daac_converter import make2dVar, makeVectorVar, makeDims, make3dVar, makeGroups, makeGlobalAttr
from emit_utils.file_checks import netcdf_ext, envi_header
from spectral.io import envi
import os
import logging


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L1B PGE outputs to DAAC compatable formats, with supporting metadata''', add_help=True)

    parser.add_argument('rdn_file', type=str, help="EMIT L1B radiance ENVI file")
    parser.add_argument('obs_file', type=str, help="EMIT L1B observation data ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('glt_file', type=str, help="EMIT L1B glt ENVI file")
    parser.add_argument('output_base', type=str, help="Output base for radiance files")
    parser.add_argument('start_datetime', type=str, help="datetime string of acquisition start in YYYMMDDHHMMSS")
    parser.add_argument('version', type=str, help="version number of processing")
    parser.add_argument('--log_file', type=str, dfault=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.log_level, filename=args.log_file)

    rdn_ds = envi.open(envi_header(args.rdn_file))
    obs_ds = envi.open(envi_header(args.obs_file))
    loc_ds = envi.open(envi_header(args.loc_file))
    glt_ds = envi.open(envi_header(args.glt_file))

    # make the netCDF4 file
    logging.info(f'Creating netCDF4 file: {netcdf_ext(args.output_base)}')
    nc_ds = Dataset(netcdf_ext(args.output_base), 'w', clobber=True, format='NETCDF4')

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

    logging.debug('Creating dimensions')
    makeDims(nc_ds, args.rdn_file, args.glt_file)

    logging.debug('Creating groups')
    makeGroups(nc_ds, ['radiance_data'])

    ## Variables (bands \ actual data)
    # instrument stuff
    logging.debug('Creating and writing sensor_band_parameters group')
    wvs = [float(d) for d in dsL1.hdr['rdn_clip']['wavelength'][0].split(',')]
    makeVectorVar(nc_ds, '/sensor_band_parameters', '/wavelength', "f4", 'number_of_bands', 'Band center wavelengths',
                  'nm', wvs)

    fwhm = [float(d) for d in dsL1.hdr['rdn_clip']['fwhm'][0].split(',')]
    makeVectorVar(nc_ds, '/sensor_band_parameters', '/fwhm', "f4", 'number_of_bands', 'full width at half max', 'nm',
                  fwhm)

    # Scan time (need to have a deeper look)
    print('Creating and writing scan_line_attributes')
    # this grabs the time stamp from the observation file
    # pdb.set_trace()
    makeVectorVar(nc_ds, '/scan_line_attributes', '/scan_start_time', "f8", 'number_of_scans', 'Scan start time (UTC)',
                  'seconds', dsL1.data['scan_start_time'])

    # GLT file
    ### In principle we can add attributes to the group instead of the variables.. e.g., map info
    # Maybe we'll do that later. for now I'm sticking to the dirty fix solution
    print('Creating and writing GLT Lookup')
    # pdb.set_trace()
    minfo = dsL1.hdr['rdn_glt']['map info']
    make2dVar(nc_ds, '/glt_data', '/Sample', "i4", 'glt_y', 'glt_x', 'GLT sample Lookup', \
              'pixel location', dsL1.data['glt sample lookup'], minfo)
    make2dVar(nc_ds, '/glt_data', '/Line', "i4", 'glt_y', 'glt_x', 'GLT line Lookup', \
              'pixel location', dsL1.data['glt line lookup'], minfo)

    # loc file
    print('Creating and writing loc_data')
    make2dVar(nc_ds, '/loc_data', '/longitude', "d", 'number_of_scans', 'pixels_per_scan',
              'Longitudes of pixel locations', 'degrees_east', dsL1.data['longitude (wgs-84)'])
    make2dVar(nc_ds, '/loc_data', '/latitude', "d", 'number_of_scans', 'pixels_per_scan',
              'Latitudes of pixel locations', 'degrees_north', dsL1.data['latitude (wgs-84)'])
    make2dVar(nc_ds, '/loc_data', '/Elevation', "d", 'number_of_scans', 'pixels_per_scan',
              'Terrain elevation at pixel locations', 'meters', dsL1.data['elevation (m)'])
    # OBS file
    print('Creating and writing obs_data')
    make2dVar(nc_ds, '/obs_data', '/path_length', "d", 'number_of_scans', 'pixels_per_scan',
              'Distance between pixel and sensor', 'meters', dsL1.data['path length (m)'])
    make2dVar(nc_ds, '/obs_data', '/sensor_zenith', "d", 'number_of_scans', 'pixels_per_scan', 'To-sensor zenith',
              'degrees', dsL1.data['to-sensor zenith (0 to 90 degrees from zenith)'])
    make2dVar(nc_ds, '/obs_data', '/to-sensor_azimuth', "d", 'number_of_scans', 'pixels_per_scan', 'To-sensor azimuth',
              'degrees', dsL1.data['to-sensor azimuth (0 to 360 degrees cw from n)'])
    make2dVar(nc_ds, '/obs_data', '/solar_zenith', "d", 'number_of_scans', 'pixels_per_scan',
              'Solar zenith angle at pixel locations', 'degrees',
              dsL1.data['to-sun zenith (0 to 90 degrees from zenith)'])
    make2dVar(nc_ds, '/obs_data', '/to-sun_azimuth', "d", 'number_of_scans', 'pixels_per_scan',
              'Solar azimuth angle at pixel locations', 'degrees',
              dsL1.data['to-sun azimuth (0 to 360 degrees cw from n)'])
    make2dVar(nc_ds, '/obs_data', '/solar_phase', "d", 'number_of_scans', 'pixels_per_scan', 'Solar phase angle',
              'degrees', dsL1.data['solar phase'])
    make2dVar(nc_ds, '/obs_data', '/slope', "d", 'number_of_scans', 'pixels_per_scan', 'Slope of the surface',
              'degrees', dsL1.data['slope'])
    make2dVar(nc_ds, '/obs_data', '/aspect', "d", 'number_of_scans', 'pixels_per_scan', 'Aspect of the surface',
              'degrees', dsL1.data['aspect'])
    make2dVar(nc_ds, '/obs_data', '/cosine_i', "d", 'number_of_scans', 'pixels_per_scan', 'cosine of SZA', 'unitless',
              dsL1.data['cosine(i)'])
    make2dVar(nc_ds, '/obs_data', '/UTC_Time', "d", 'number_of_scans', 'pixels_per_scan', 'Slope of the surface',
              'Decimal Hours', dsL1.data['utc time'])
    make2dVar(nc_ds, '/obs_data', '/Earth-sun-distance', "d", 'number_of_scans', 'pixels_per_scan',
              'Distance between Earth and Sun', 'Astronomical Distance Units', dsL1.data['earth-sun distance (au)'])
    # Radiance file
    print('Creating and writing rdn_data')
    make3dVar(nc_ds, '/rdn_data', '/rdn', "f4", 'number_of_bands', 'number_of_scans', 'pixels_per_scan',
              'At-sensor Radiance', 'uW/cm/SR/nm', dsL1.data['rdn'])

    # Finished
    nc_ds.close()
    print('Successfully created', nc_fullpath)

    return


