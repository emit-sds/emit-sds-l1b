import argparse
import os
import uuid

import cv2 as cv
import numpy as np
from osgeo import gdal
import spectral
from spectral.io import envi

def single_image_ortho(img_dat, in_glt, glt_nodata_value=0):
    glt = in_glt.copy()
    outdat = np.zeros((glt.shape[0], glt.shape[1], img_dat.shape[-1])) - 9999
    valid_glt = np.all(glt != glt_nodata_value, axis=-1)
    glt[valid_glt] -= 1
    outdat[valid_glt, :] = img_dat[glt[valid_glt, 1], glt[valid_glt, 0], :]
    return outdat

def create_mosaic(input_files, output_file):
    vrt_options = gdal.BuildVRTOptions(resolution='average')
    vrt_ds = gdal.BuildVRT(f'/vsimem/{uuid.uuid4().hex}.vrt', input_files, options=vrt_options)
    translate_options = gdal.TranslateOptions(format="COG", creationOptions=["TILING_SCHEME=GoogleMapsCompatible"])
    gdal.Translate(output_file, vrt_ds, options=translate_options)

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

    waves = [660, 550, 440]
    rgb = []
    indices = {}
    start = 0
    for file in args.rdn:
        img = spectral.envi.open(file.replace('.img', '.hdr'))
        indices[os.path.basename(file)] = [start, start + img.nrows]
        start += img.nrows
        wavelengths = np.array(img.metadata.get('wavelength', []), dtype=float)
        idx = [np.argmin(np.abs(wavelengths - wave)) for wave in waves]
        rgb.append(img.read_bands(idx))

    rgb = np.vstack(rgb)
    mask = rgb[:, :, 0] == args.nodata
    rgb[mask] = np.nan
    rgb_mask = rgb.copy()
    rgb_mask[rgb_mask > 50] = np.nan
    rgb = np.clip(rgb, 0, 50)
    lower = np.nanmin(rgb_mask, axis=(0, 1))
    upper = np.nanmax(rgb_mask, axis=(0, 1))
    rgb = (rgb - lower) / (upper - lower) * 254
    rgb[mask] = 0
    rgb[np.isnan(rgb)] = 0
    rgb = rgb.astype('uint8')

    clip = 10.0
    tile = (16, 16)
    rgb_adj = rgb.copy()
    clahe = cv.createCLAHE(clipLimit=clip, tileGridSize=tile)
    for b in range(rgb_adj.shape[-1]):
        rgb_adj[:, :, b] = clahe.apply(rgb_adj[:, :, b])

    ort_rgbfiles = []
    for rdn_file, glt_file in zip(args.rdn, args.glt):
        glt_dataset = envi.open(glt_file.replace('.img', '.hdr'))
        glt = glt_dataset.open_memmap(writeable=False, interleave='bip').copy()
        del glt_dataset
        glt_dataset = gdal.Open(glt_file)
        start, end = indices[os.path.basename(rdn_file)]
        rgb_ort = single_image_ortho(rgb_adj[start:end], glt)
        rgb_ort[rgb_ort == args.nodata] = 0
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
        ort_rgbfiles.append(rgb_dset)

    create_mosaic(ort_rgbfiles, args.output)

if __name__ == '__main__':
    main()
