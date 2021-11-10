# David R Thompson

import numpy as np
import pylab as plt
from glob import glob
from spectral.io import envi
from scipy.linalg import solve, inv

parser = argparse.ArgumentParser(description='Remove Dyson Ghost.')
parser.add_argument('input', type=str, help='input file')
parser.add_argument('ghost_matrix', type=str, help='ghost matrix')
parser.add_argument('--center', type=float, default=650.0, 
   help='center of symmetry')
parser.add_argument('output', type=str, help='output file')
args = parser.parse_args()

inimg = envi.open(args.input+'.hdr')
meta = inimg.metadata.copy()
frames = inimg.load()
ghost = np.squeeze(envi.open(args.ghost_matrix+'.hdr').load())
left, right, short, long = 25, 1265, 21, 314
nrow, ncol = 480, 1280
center = args.center

A = np.concatenate((np.concatenate((np.eye(nrow),ghost.T),axis=1),
                    np.concatenate((ghost.T,np.eye(nrow)),axis=1)), axis=0)
Ainv = inv(A)
out = []

for frame_ind in range(frames.shape[0]):
    frame = np.squeeze(frames[frame_ind,:,:]).T
    new = np.zeros(frame.shape)
    for col in range(15,int(center)):
        tcol = int(np.ceil(center+(center-col))) # APPROXIMATION
        b = np.concatenate((frame[:,col],frame[:,tcol]),axis=0)
        x = Ainv @ b[:,np.newaxis]
        new[:,col] = x[:nrow,0]
        new[:,tcol] = x[nrow:,0]
    out.append(new)

out = np.array(out,dtype=np.float32)
envi.save_image('.hdr',out,metadata=meta,ext='',force=True)
print('done')

