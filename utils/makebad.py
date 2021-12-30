# David R Thompson
import argparse, sys, os
import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.stats import norm
from scipy.linalg import solve, inv
from astropy import modeling
from sklearn.linear_model import RANSACRegressor
from skimage.filters import threshold_otsu
import json


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

# bad row, column ranges (inclusive, zero indexed)
manual_bads=[(1271,(5,6)),
  (1258,198),
  (1256,193),
  (1231,63),
  (1228,180),
  (1194,205),
  (1192,(86,87)),
  (1192,63),
  (1177,125),
  (1176,125),
  (1173,(130,131)),
  (1172,(130,131)),
  ((1168,1170),(40,43)),
  (1165,68),
  (1164,310),
  (1164,383),
  (1243,211),
  (1141,(315,316)),
  (1141,73),
  (1140,312),
  (1139,312),
  (1138,(311,313)),
  (1123,213),
  (1123,75),
  (1109,78),
  (1108,78),
  (1107,232),
  (1077,(217,218)),
  (1076,(217,218)),
  (1067,129),
  (1059,148),
  (1033,113),
  (1009,231),
  (1005,114),
  (1001,266),
  (1000,49),
  (999,(93,95)),
  (998,(93,95)),
  (988,104),
  (987,107),
  (963,280),
  (963,281),
  (963,34),
  (962,(281,283)),
  (955,(103,105)),
  (954,(103,105)),
  (951,270),
  (951,177),
  (951,268),
  (950,(268,269)),
  (948,126),
  (948,270),
  (947,269),
  (946,106),
  (945,343),
  (944,173),
  (940,59),
  (923,88),
  (922,252),
  (913,166),
  (899,105),
  (880,64),
  (858,424),
  ((828,844),(436,455)),
  (828,268),
  (827,213),
  (821,226),
  (815,186),
  (809,79),
  (805,(134,136)),
  (803,(80,81)),
  (803,(134,136)),
  (802,(80,81)),
  ((769,781),(409,416)),
  (794,(362,363)),
  (763,49),
  (762,(49,50)),
  (761,49),
  (761,122),
  (760,(121,122)),
  (752,158),
  (752,163),
  (746,448),
  (729,311),
  (721,(247,248)),
  (719,249),
  (719,245),
  (678,401),
  (645,215),
  (644,215),
  (633,324),
  (625,53),
  (606,31),
  (591,416),
  (584,416),
  (576,(144,146)),
  (575,(144,146)),
  (574,(144,146)),
  (573,353),
  (569,364),
  (568,364),
  (532,83),
  (529,350),
  (507,(253,254)),
  (506,(253,254)),
  (505,(253,254)),
  (490,257),
  (485,100),
  (462,327),
  (461,78),
  (447,265),
  (446,(264,266)),
  (445,265),
  (442,217),
  (434,(201,202)),
  (425,258),
  (424,(257,259)),
  (423,258),
  (376,129),
  (374,45),
  (363,119),
  (356,242),
  (349,225),
  (348,(34,35)),
  (346,219), 
  (345,218),
  (345,(218,220)), 
  (344,(217,219)), 
  (343,218),
  (325,(168,169)), 
  (310,237),
  (308,237),
  (307,468),
  (301,440),
  (295,148),
  (238,54),
  (233,102),
  (233,(62,63)),
  (232,102),
  (231,(100,102)),
  (227,(99,100)),
  (226,(99,100)),
  (225,99),
  (217,241),
  (216,(240,242)),
  (215,(236,242)),
  (214,(236,237)),
  (213,(236,238)),
  (212,(244,246)),
  (212,238),
  (211,(244,246)),
  (211,239),
  (211,(241,243)),
  (210,(240,245)),
  (206,336),
  (202,401),
  (201,461),
  (196,251),
  (195,251),
  (194,247),
  (193,247),
  (191,235),
  (189,110),
  (188,237),
  (184,130),
  (179,243),
  (167,247),
  (157,188),
  (136,133),
  (135,133),
  (125,154),
  (119,54),
  (118,249),
  (117,157),
  (117,87),
  (115,41),
  (108,355),
  (108,255),
  (104,135),
  (104,104),
  (89,479),
  (87,176),
  (86,(175,177)),
  (85,176),
  (81,390),
  (72,126),
  (54,(224,225)),
  (53,(224,225)),
  (52,255),
  (52,154),
  (52,(254,255)),
  (51,255),
  (48,109),
  (47,153),
  (46,153),
  (26,273),
  (27,26),
  (3,114)]

           

def main():

    description = "Calculate Flat field"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--cue_channel',default=148,type=int)
    parser.add_argument('--ref_lo',default=99,type=int)
    parser.add_argument('--ref_hi',default=1180,type=int)
    parser.add_argument('--hw_lo',default=50,type=int)
    parser.add_argument('--hw_hi',default=180,type=int)
    parser.add_argument('--selection',type=str,default='spatial')
    parser.add_argument('--badmap_out',type=str,default=None)
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))
 
    if int(infile.metadata['data type']) == 2:
        dtype = np.uint16
    elif int(infile.metadata['data type']) == 4:
        dtype = np.float32
    else:
        raise ValueError('Unsupported data type')
    if infile.metadata['interleave'] != 'bil':
        raise ValueError('Unsupported interleave')


    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    lines = int(infile.metadata['lines'])
    nframe = rows * columns
    margin=2

    flat  = np.zeros((rows,columns))
    count = np.zeros((rows,columns))
    sumsq = np.zeros((rows,columns))
    ref = np.zeros((lines,columns))
    allctrs,alllines = [],[]
    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            ref[line,:] = frame[args.cue_channel, :]

    thresh = np.sort(ref,axis=0)
    thresh = thresh[-10,:]

    with open(args.input,'rb') as fin:

        for line in range(lines):

            # Read a frame of data
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            reference = frame[args.cue_channel, :]
            use = np.where(reference>thresh)[0]

            print(line,len(use),np.median(use),thresh)
            flat[:,use] = flat[:,use] + frame[:,use] 
            count[:,use] = count[:,use] + 1
            sumsq[:,use] = sumsq[:,use] + pow(frame[:,use],2)

        mean_sumsq = sumsq / count
        flat = flat / count

        rowmean = flat[:,30:1250].mean(axis=1)
        rowstdev = flat[:,30:1250].std(axis=1)
        stdev = np.sqrt(mean_sumsq - pow(flat,2))
        stdev[np.logical_not(np.isfinite(stdev))] = 0
        bad = np.logical_or(np.logical_or(stdev==0,
              (abs(flat.T-rowmean)>rowstdev*20).T),stdev>100)
    
    bad[:,:25] = 0
    bad[:,1265:] = 0

    for bad_cols, bad_rows in manual_bads:
        if type(bad_rows)==int:
            rows_range = [bad_rows]
        else:
            rows_range = range(bad_rows[0],bad_rows[1]+1)
        if type(bad_cols)==int:
            cols_range = [bad_cols]
        else:
            cols_range = range(bad_cols[0],bad_cols[1]+1)
        for col in cols_range:
            for row in rows_range:
                bad[row,col] = 1
   #plt.hist(stdev.flatten(),500)
   #plt.figure()
   #plt.imshow(bad)
   #plt.show()
    bads = 0
    bad_map = bad.copy()
    bad_map = np.array(bad_map,dtype=np.int16)
    for column in range(bad_map.shape[1]):
        state_machine = 0
        for row in range(bad_map.shape[0]):
            if bad[row,column]:
                state_machine = state_machine + 1
                bad_map[row,column] = -state_machine
                print(row,column,state_machine)
                bads = bads + 1
            else:
                state_machine = 0
    print('total bads:',bads)
    bad_map = bad_map.reshape((rows,columns,1))
    envi.save_image(args.output+'.hdr',
        bad_map, interleave='bsq', ext='', force=True)

if __name__ == '__main__':

    main()
