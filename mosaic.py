import argparse
import os
import uuid

import cv2 as cv
import numpy as np
from osgeo import gdal
import spectral
from spectral.io import envi

def process_glt(glt_file, ulx, lrx):
    ds = gdal.Open(glt_file)
    band1 = ds.GetRasterBand(1).ReadAsArray()
    gt = ds.GetGeoTransform()
    gp = ds.GetProjection()

    width = ds.RasterXSize
    height = ds.RasterYSize
    uly = gt[3]
    lry = uly + gt[5] * height
    lon = np.linspace(ulx + gt[1] / 2, lrx - gt[1] / 2, width)
    lat = np.linspace(uly + gt[5] / 2, lry - gt[5] / 2, height)
    return gt, gp, band1, lon, lat

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

def array_to_gdal(rgb_ort, proj, geotrans, nodata = -9999):
    rgb_ort[rgb_ort == nodata] = 0
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
    mask = rgb[:, :, 0] == -9999
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
        rgb_adj[mask.sum(axis=1)==0,:, b] = clahe.apply(rgb_adj[mask.sum(axis=1)==0,:, b])
        zero_but_data  = (rgb_adj[...,b] == 0) * ~mask
        rgb_adj[zero_but_data,b] += 1

    west_files = []
    east_files = []

    for rdn_file, glt_file in zip(args.rdn, args.glt):

        info = gdal.Info(glt_file, format='json')
        ulx = info['cornerCoordinates']['upperLeft'][0]
        lrx = info['cornerCoordinates']['lowerRight'][0]
        print(f'File: {glt_file}')
        print(f'  Upper Left Longitude: {ulx}')
        print(f'  Lower Right Longitude: {lrx}')

        geotrans, geoproj, band1, lon, lat = process_glt(glt_file, ulx, lrx)
        glt_dataset = envi.open(glt_file.replace('.img', '.hdr'))
        glt = glt_dataset.open_memmap(writeable=False, interleave='bip').copy()
        del glt_dataset

        start, end = indices[os.path.basename(rdn_file)]

        if ulx < 180 and lrx > 180:
            print(f'{os.path.basename(rdn_file)} Crosses antimeridian, splitting\n')

            lright = np.argwhere(lon >= 179.999)[0][0]
            lband = band1[:, :lright + 1]
            ltop, lbottom = np.argwhere(np.sum(lband == 0, axis=1) != lband.shape[1])[[0, -1]].flatten()
            glt_left = glt[ltop:lbottom + 1, 0:lright + 1]

            rleft = np.argwhere(lon >= 180.001)[0][0]
            rright = band1.shape[1] - 1
            rband = band1[:, rleft:]
            rtop, rbottom = np.argwhere(np.sum(rband == 0, axis=1) != rband.shape[1])[[0, -1]].flatten()
            glt_right = glt[rtop:rbottom + 1, rleft:rright + 1]
            gt_right = list(geotrans)
            gt_right[0] = -180 - (180 - lon[rleft]) if lon[rleft] > 180 else lon[rleft]
            gt_right[3] = lat[rtop]

            rgb_ort_left = single_image_ortho(rgb_adj[start:end], glt_left).astype(int)
            rgb_ort_right = single_image_ortho(rgb_adj[start:end], glt_right).astype(int)

            east_files.append(array_to_gdal(rgb_ort_left,
                                          geoproj,
                                          geotrans))

            west_files.append(array_to_gdal(rgb_ort_right,
                                          geoproj,
                                          gt_right))

        else:

            rgb_ort = single_image_ortho(rgb_adj[start:end], glt).astype(int)

            if ulx > 180:
                print(f'{os.path.basename(rdn_file)}  Lies east of east antimeridian, moving west\n')
                geotrans[0] -= 360

            gdal_dset = array_to_gdal(rgb_ort,
                            geoproj,
                            geotrans)

            if ulx > 180 or lrx <= 0:
                print(f'{os.path.basename(rdn_file)}  Assigned to: WEST\n')
                west_files.append(gdal_dset)
            else:
                print(f'{os.path.basename(rdn_file)}  Assigned to: EAST\n')
                east_files.append(gdal_dset)


    crosses_antimeridian = bool(west_files and east_files and west_files != east_files)
    print(f'Crosses antimeridian: {crosses_antimeridian}\n')

    if crosses_antimeridian:
        output_base, ext = os.path.splitext(args.output)
        east_output = f'{output_base}_2{ext}'

        print(f'Creating separate mosaics:\n  West: {args.output}\n  East: {east_output}\n')

        if west_files:
            print(f'Building west mosaic from {len(west_files)} files')
            create_mosaic(west_files, args.output)

        if east_files:
            print(f'Building east mosaic from {len(east_files)} files')
            create_mosaic(east_files, east_output)
    else:
        print('Creating single mosaic (no antimeridian crossing or only one input)\n')
        create_mosaic(west_files + east_files, args.output)

if __name__ == '__main__':
    main()
