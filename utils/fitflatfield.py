#!/usr/bin/env python
# David R Thompson 
# Spring 2015
# Jet Propulsion Laboratory, California Institute of Technology

import os, sys, argparse, time
import spectral
from scipy.optimize import minimize
import numpy as np
import spectral.io.envi as envi
from numba import jit
from emit2dark import bad_flag, dark_from_file

def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


# Squared error with an L2-norm regularizer.  This operates on just one
# spectral channel.  It requires the images' sufficient statistics in that
# channel to be precomputed and stored in "selfsq" and "adj."  The error
# takes the form sum((I(r,c+1) - I(r,c))^2), e.g. squared difference of 
# horizontally-adjacent pixels. 
#
#   ff - a vector of size [nsamples], 1 gain value per cross-track element
#   selfsq - a vector of size [nsamples], the column-sums of squared radiances
#   self - a vector of size [nsamples], the column-sums of radiances
#   adj - a vector of size [nsamples-1], the column-sums of products of 
#         horizontally-adjacent image elements. 
#   reg - a regularizer value, typically 0-1
#@jit
def err(x, self, selfsq, adj, reg, scalesq):
  nc = len(self)
  B = x[:nc] # dark level (additive)
  A = x[nc:] # flat field (multiplicative)
  er = 0
  for c in range(nc-1):
    Q1  = selfsq[c]*A[c]*A[c]
    Q2  = self[c]*A[c]*B[c] 
    Q3  = adj[c]*A[c]*A[c+1]
    Q4  = self[c]*A[c]*B[c+1] 
    Q5  = B[c]*B[c]
    Q6  = self[c+1]*A[c+1]*B[c]
    Q7  = B[c+1]*B[c]
    Q8  = selfsq[c+1]*A[c+1]*A[c+1]
    Q9  = self[c+1]*A[c+1]*B[c+1]
    Q10 = B[c+1]*B[c+1]
    er = er + Q1+2*Q2-2*Q3-2*Q4+Q5-2*Q6-2*Q7+Q8+2*Q9+Q10
  nz = float(nc*2-1.0)
  er = er/nz/scalesq + reg[0]*(pow(A-1.0,2).sum()) + reg[1]*(pow(B,2).sum())
  return er


# Gradient of the error has a pleasant closed form
#@jit
def grad(x, self, selfsq, adj, reg, scalesq):

  nc = len(self)
  B = x[:nc] # dark level (additive)
  A = x[nc:] # flat field (multiplicative)
  partialA = np.zeros(A.shape)
  partialB = np.zeros(B.shape)
  for c in range(nc-1):
    Q1_A  = selfsq[c]*2*A[c]
    Q1_Ap = 0
    Q1_B  = 0
    Q1_Bp = 0
    Q2_A  = self[c]*B[c] 
    Q2_Ap = 0
    Q2_B  = self[c]*A[c] 
    Q2_Bp = 0
    Q3_A  = adj[c]*A[c+1]
    Q3_Ap = adj[c]*A[c]
    Q3_B  = 0
    Q3_Bp = 0
    Q4_A  = self[c]*B[c+1] 
    Q4_Ap = 0
    Q4_B  = 0
    Q4_Bp = self[c]*A[c] 
    Q5_A  = 0
    Q5_Ap = 0
    Q5_B  = 2*B[c]
    Q5_Bp = 0
    Q6_A  = 0
    Q6_Ap = self[c+1]*B[c]
    Q6_B  = self[c+1]*A[c+1]
    Q6_Bp = 0
    Q7_A  = 0
    Q7_Ap = 0
    Q7_B  = B[c+1]
    Q7_Bp = B[c]
    Q8_A  = 0
    Q8_Ap = selfsq[c+1]*2*A[c+1]
    Q8_B  = 0
    Q8_Bp = 0
    Q9_A  = 0
    Q9_Ap = self[c+1]*B[c+1]
    Q9_B  = 0
    Q9_Bp = self[c+1]*A[c+1]
    Q10_A  = 0
    Q10_Ap = 0
    Q10_B  = 0
    Q10_Bp = 2*B[c+1]
    partialA[c] = partialA[c]+Q1_A+2*Q2_A-2*Q3_A-2*Q4_A+Q5_A-2*Q6_A-\
                                2*Q7_A+Q8_A+2*Q9_A+Q10_A
    partialB[c] = partialB[c]+Q1_B+2*Q2_B-2*Q3_B-2*Q4_B+Q5_B-2*Q6_B-\
                                2*Q7_B+Q8_B+2*Q9_B+Q10_B
    partialA[c+1] = partialA[c+1]+Q1_Ap+2*Q2_Ap-2*Q3_Ap-2*Q4_Ap+Q5_Ap-2*Q6_Ap-\
                                2*Q7_Ap+Q8_Ap+2*Q9_Ap+Q10_Ap
    partialB[c+1] = partialB[c+1]+Q1_Bp+2*Q2_Bp-2*Q3_Bp-2*Q4_Bp+Q5_Bp-2*Q6_Bp-\
                                2*Q7_Bp+Q8_Bp+2*Q9_Bp+Q10_Bp

  # regularizer
  nz = float(nc*2-1.0)
  partialA = partialA/nz/scalesq + reg[0]*2*(A-1.0)
  partialB = partialB/nz/scalesq + reg[1]*2*(B)
  
  # pack additive and multiplicative terms into a single state vector
  partialX = np.zeros(x.shape)
  partialX[:nc] = partialB
  partialX[nc:] = partialA
  return partialX


def main():

    description = 'Correct striping in pushbroom radiance data';
    parser = argparse.ArgumentParser()

    # Required 
    parser.add_argument('input', help='Input radiance image')
    parser.add_argument('output', help='Output radiance image')
    parser.add_argument('flatfield', help='New (secondary) flat field')
    parser.add_argument('darkfield', help='New (secondary) dark field')
    parser.add_argument('--reg','-r',help='Regularizer', default=0.0001) 
    parser.add_argument('--max_radiance',help='Max. Radiance', default=10) 
    parser.add_argument('--margin',help='Spatial margin', default=50) 
    args = vars(parser.parse_args())

    # Define local variables
    inhdr  = find_header(args['input'])
    outhdr = args['output'] + '.hdr'
    ffhdr  = args['flatfield']+'.hdr'
    dkhdr  = args['darkfield']+'.hdr'
    reg = np.array([float(args['reg']), float(args['reg'])])

    # First I/O pass: accumulate "sufficient statistics" for FF optimization.
    I = envi.open(inhdr)
    nrows, ncols, nbands = I.nrows, I.ncols, I.nbands
    meta = I.metadata.copy()
    self   = np.zeros((nbands,ncols))
    selfsq = np.zeros((nbands,ncols))
    adj    = np.zeros((nbands,ncols-1))
    scaling   = None
    nused = np.zeros((nbands,ncols))

    with open(args['input'],'rb') as fin:

        frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)

        while frame.size>0:
            frame = frame.reshape((nbands, ncols)) # BIL
            if any((frame<-9989).flatten()):
               continue

            notbright = np.max(frame, axis=0) < args['max_radiance']
            valid = np.all(frame > bad_flag, axis=0)
            use = np.logical_and(notbright, valid)
            
            # Cached values are stored in arrays the size of the active FPA region
            # We sum over all lines (and will later turn this into a mean)
            # The final values will be based on the mean squared error at each
            # cross-track discontinuity
            if sum(use)>1: 
                 self[:,use]   = self[:,use] + frame[:,use]
                 selfsq[:,use] = selfsq[:,use] + pow(frame[:,use],2)
                 diffs = (frame[:,:-1] * frame[:,1:])
                 adj[:,use[:-1]]    = adj[:,use[:-1]] + diffs[:,use[:-1]] 
                 nused[:,use] = nused[:,use]+1
            frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)

    self   = self   / nused
    selfsq = selfsq / nused
    adj    = adj    / nused[:,:-1]
    scaling = self[:, args.margin:-args.margin].mean(axis=1)

    # Second pass: optimize and correct
    ff = np.ones((nbands, ncols))
    dk = np.zeros((nbands, ncols))
    nc = self.shape[1]
    
    # Optimize each band independently.  Unexpected out-of-range values can 
    # cause convergence failure, but 100 iterations is enough for most cases
    for b in range(nbands): 
      print('Band %i'%(b));  
      f0 = np.r_[np.zeros((nc,1)), np.ones((nc,1))]
      arg = (self[b,:], selfsq[b,:], adj[b,:], reg, np.array([scaling[b]*scaling[b]]))
      opts = {'maxiter':1000, 'disp':True}
      res = minimize(err, f0, method='CG',args=arg, jac=grad, options=opts, tol=1e-6)
      if res.success:
          dk[b,:] = res.x[:nc]
          ff[b,:] = res.x[nc:]
      else:
        print(res.message)

    # Write the full-sized optimized and mean-divided flat field
    ff3d = np.array(ff.reshape((nbands,ncols,1)), dtype=np.float32)
    envi.save_image(ffhdr, ff3d, force=True, ext="", interleave='bsq')

    # Write the full-sized optimized and mean-divided dark field
    dk3d = np.array(dk.reshape((nbands, ncols,1)), dtype=np.float32)
    envi.save_image(dkhdr, dk3d, force=True, ext="", interleave='bsq')

    # Second I/O pass applies the correction to the whole image,
    # one frame at a time
    Icorr = envi.create_image(outhdr, meta, force=True, ext="")

    with open(args['input'],'rb') as fin:
        with open(args['output'],'wb') as fout:

            frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)

            while frame.size > 0:

                 frame = frame.reshape((nbands, ncols)) # BIL
                 if frame.shape[0]<1:
                     continue
                 new_frame = dk+frame*ff # It's that easy.
                 new_frame = np.array(new_frame, dtype=np.float32)  
                 new_frame.tofile(fout)
                 frame = np.fromfile(fin, count=nbands*ncols, dtype=np.float32)
                   
if __name__ == '__main__':
    main()


