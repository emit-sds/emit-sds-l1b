#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import spectral
from spectral.io import envi
import uuid
from osgeo import gdal


def create_mosaic(input_files, output_file):
    vrt_options = gdal.BuildVRTOptions(resolution='average')
    vrt_ds = gdal.BuildVRT(f'/vsimem/{uuid.uuid4().hex}.vrt', input_files, options=vrt_options)
    translate_options = gdal.TranslateOptions(format="COG",
                                              creationOptions=["TILING_SCHEME=GoogleMapsCompatible"])
    gdal.Translate(output_file, vrt_ds, options=translate_options)


def stretch_statistics(input_files, nodata=-9999):
    waves = [660, 550, 440]
    lower = {'r': [], 'g': [], 'b': []}
    upper = {'r': [], 'g': [], 'b': []}
    for file in input_files:
        img = spectral.envi.open(file.replace('.img', '.hdr'), file)
        wavelengths = np.array(img.metadata.get('wavelength', []), dtype=float)
        cube = img.load()
        cube = np.where(cube == nodata, np.nan, cube)
        for i, wave in enumerate(waves):
            idx = np.argmin(np.abs(wavelengths - wave))
            band = cube[:, :, idx]
            lower['rgb'[i]].append(np.nanpercentile(band, 1))
            upper['rgb'[i]].append(np.nanpercentile(band, 99))
    stats = {'lower_rgb': [], 'upper_rgb': []}
    for band in 'rgb':
        stats['lower_rgb'].append(float(np.nanmedian(lower[band])))
        stats['upper_rgb'].append(float(np.nanmedian(upper[band])))
    return stats

# TODO import from apply_glt in emit-utils
def single_image_ortho(img_dat, in_glt, glt_nodata_value=0):
    glt = in_glt.copy()
    outdat = np.zeros((glt.shape[0], glt.shape[1], img_dat.shape[-1])) - 9999
    valid_glt = np.all(glt != glt_nodata_value, axis=-1)
    glt[valid_glt] -= 1
    outdat[valid_glt, :] = img_dat[glt[valid_glt, 1], glt[valid_glt, 0], :]
    return outdat


def create_quicklook(input_file, glt_file, stats):
    img = spectral.envi.open(input_file.replace('.img', '.hdr'), input_file)
    wavelengths = np.array(img.metadata.get('wavelength', []), dtype=float)
    cube = img.load()
    rgb_waves = [660, 550, 440]
    rgb = []
    for i, wave in enumerate(rgb_waves):
        idx = np.argmin(np.abs(wavelengths - wave))
        band = cube[:, :, idx].astype(float)
        mask = band == -9999
        band = np.clip((band - stats['lower_rgb'][i]) / (stats['upper_rgb'][i] - stats['lower_rgb'][i]), 0, 1)
        rgb.append(band)
    rgb = (np.dstack(rgb) * 254).astype(np.uint8)
    rgb[~mask[:, :, 0]] += 1
    glt_dataset = envi.open(glt_file.replace('.img', '.hdr'))
    glt = glt_dataset.open_memmap(writeable=False, interleave='bip').copy()
    del glt_dataset
    glt_dataset = gdal.Open(glt_file)
    rgb_ort = single_image_ortho(rgb, glt)
    rgb_ort[rgb_ort == -9999] = 0
    proj = glt_dataset.GetProjection()
    geotrans = glt_dataset.GetGeoTransform()
    driver = gdal.GetDriverByName('MEM')
    rows, cols, bands = rgb_ort.shape
    rgb_dset = driver.Create('', cols, rows, bands, gdal.GDT_Byte)
    rgb_dset.SetProjection(proj)
    rgb_dset.SetGeoTransform(geotrans)
    for i in range(bands):
        band = rgb_dset.GetRasterBand(i + 1)
        band.WriteArray(rgb_ort[:, :, i].astype(np.uint8))
        band.SetNoDataValue(0)
    return rgb_dset


def main():
    parser = argparse.ArgumentParser(description="Create RGB mosaic from L1B radiance")
    parser.add_argument('--rdn', nargs='+', required=True, help='L1B radiance files')
    parser.add_argument('--glt', nargs='+', required=True, help='GLT files')
    parser.add_argument('--output', required=True, help='Output mosaic file')
    parser.add_argument('--nodata', type=float, default=-9999, help='Nodata value')
    args = parser.parse_args()

    if len(args.rdn) != len(args.glt):
        raise ValueError("The number of radiance files and GLT files must be equal")

    args.rdn.sort()
    args.glt.sort()
    stats = stretch_statistics(args.input, nodata=args.nodata)
    ort_rgbfiles = []
    for input_file, glt_file in zip(args.rdn, args.glt):
        ort_rgbfiles.append(create_quicklook(input_file, glt_file, stats))
    create_mosaic(ort_rgbfiles, args.output)


if __name__ == '__main__':
    main()
