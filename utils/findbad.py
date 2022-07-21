# David R Thompson
from __future__ import division, print_function
import numpy as np
from sklearn.covariance import MinCovDet
from spectral.io import envi
import os, sys, argparse
from fpa import FPA
from scipy.linalg import inv
import pylab as plt
import logging
from numpy.random import permutation
import ray



# From https://github.com/dganguli/robust-pca
class R_pca:

    def __init__(self, D, mu=None, lmbda=None):
        self.D = D
        self.S = np.zeros(self.D.shape)
        self.Y = np.zeros(self.D.shape)

        if mu:
            self.mu = mu
        else:
            self.mu = np.prod(self.D.shape) / (4 * np.linalg.norm(self.D, ord=1))

        self.mu_inv = 1 / self.mu

        if lmbda:
            self.lmbda = lmbda
        else:
            self.lmbda = 1 / np.sqrt(np.max(self.D.shape))

    @staticmethod
    def frobenius_norm(M):
        return np.linalg.norm(M, ord='fro')

    @staticmethod
    def shrink(M, tau):
        return np.sign(M) * np.maximum((np.abs(M) - tau), np.zeros(M.shape))

    def svd_threshold(self, M, tau):
        U, S, V = np.linalg.svd(M, full_matrices=False)
        return np.dot(U, np.dot(np.diag(self.shrink(S, tau)), V))

    def fit(self, tol=None, max_iter=1000, iter_print=100):
        iter = 0
        err = np.Inf
        Sk = self.S
        Yk = self.Y
        Lk = np.zeros(self.D.shape)

        if tol:
            _tol = tol
        else:
            _tol = 1E-7 * self.frobenius_norm(self.D)

        #this loop implements the principal component pursuit (PCP) algorithm
        #located in the table on page 29 of https://arxiv.org/pdf/0912.3599.pdf
        while (err > _tol) and iter < max_iter:
            Lk = self.svd_threshold(
                self.D - Sk + self.mu_inv * Yk, self.mu_inv)                            #this line implements step 3
            Sk = self.shrink(
                self.D - Lk + (self.mu_inv * Yk), self.mu_inv * self.lmbda)             #this line implements step 4
            Yk = Yk + self.mu * (self.D - Lk - Sk)                                      #this line implements step 5
            err = self.frobenius_norm(self.D - Lk - Sk)
            iter += 1
            if (iter % iter_print) == 0 or iter == 1 or iter > max_iter or err <= _tol:
                print('iteration: {0}, error: {1}'.format(iter, err))

        self.L = Lk
        self.S = Sk
        return Lk, Sk

    def plot_fit(self, size=None, tol=0.1, axis_on=True):

        n, d = self.D.shape

        if size:
            nrows, ncols = size
        else:
            sq = np.ceil(np.sqrt(n))
            nrows = int(sq)
            ncols = int(sq)

        ymin = np.nanmin(self.D)
        ymax = np.nanmax(self.D)
        print('ymin: {0}, ymax: {1}'.format(ymin, ymax))

        numplots = np.min([n, nrows * ncols])
        plt.figure()

        for n in range(numplots):
            plt.subplot(nrows, ncols, n + 1)
            plt.ylim((ymin - tol, ymax + tol))
            plt.plot(self.L[n, :] + self.S[n, :], 'r')
            plt.plot(self.L[n, :], 'b')
            if not axis_on:
                plt.axis('off')



def conditional_gaussian(mu: np.array, C: np.array, window: np.array, remain: np.array, x: np.array) -> \
        (np.array, np.array):
    """Define the conditional Gaussian distribution for convenience.

    len(window)+len(remain)=len(x)

    Args:
        mu: mean values
        C: matrix for conditioning
        window: contains all indices not in remain
        remain: contains indices of the observed part x1
        x: values to condition with

    Returns:
        (np.array, np.array): conditional mean, conditional covariance

    """
    w = np.array(window)[:,np.newaxis]
    r = np.array(remain)[:,np.newaxis]
    C11 = C[r, r.T]
    C12 = C[r, w.T]
    C21 = C[w, r.T]
    C22 = C[w, w.T]

    Cinv = svd_inv(C11)
    conditional_mean = mu[window] + C21 @ Cinv @ (x-mu[remain])
    conditional_cov = C22 - C21 @ Cinv @ C12
    return conditional_mean, conditional_cov


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')



def conditional_gaussian(mu: np.array, C: np.array, window: np.array, 
      remain: np.array, x: np.array) -> \
        (np.array, np.array):
    """Define the conditional Gaussian distribution for convenience.

    len(window)+len(remain)=len(x)

    Args:
        mu: mean values
        C: matrix for conditioning
        window: contains all indices not in remain
        remain: contains indices of the observed part x1
        x: values to condition with

    Returns:
        (np.array, np.array): conditional mean, conditional covariance

    """
    w = np.array(window)[:,np.newaxis]
    r = np.array(remain)[:,np.newaxis]
    C11 = C[r, r.T]
    C12 = C[r, w.T]
    C21 = C[w, r.T]
    C22 = C[w, w.T]

    Cinv = inv(C11) # Could make better use of PSD matrix here
    conditional_mean = mu[window] + C21 @ Cinv @ (x-mu[remain])
    conditional_cov = C22 - C21 @ Cinv @ C12
    return conditional_mean, conditional_cov



def main():

    description = "Fix bad pixels"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--num_cpus',type=int,default=20)
    parser.add_argument('--frames',type=int,default=3)
    parser.add_argument('--thresh',type=float,default=5)
    parser.add_argument('--npca',type=float,default=100)
    parser.add_argument('config')
    parser.add_argument('output_bad')
    parser.add_argument('output')
    args = parser.parse_args()

    fpa = FPA(args.config)

    ray.init(num_cpus=args.num_cpus)

    infile = envi.open(find_header(args.input))
    print(fpa.bad_element_file)
    badfile = envi.open(find_header(fpa.bad_element_file))
    bad = np.squeeze(badfile.load()[:,:,0])

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

    envi.write_envi_header(args.output+'.hdr',infile.metadata)

    to_use = set(np.array(np.linspace(0,lines-1,args.frames),dtype=int))
    print(to_use)

    # First pass reading selected frames
    frames = []
    with open(args.input,'rb') as fin:
       for line in range(lines):
           # Read a frame of data
           if line%10==0:
               logging.info('Line '+str(line))
           frame = np.fromfile(fin, count=nframe, dtype=dtype)
           if not line in to_use:
               continue
           print(line)
           frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
           frames.append(frame.T)
       spectra = np.array(frames).reshape((columns*len(frames),rows))
    
    # Next we generate a robust mean and covariance matrix for the image
    model = R_pca(spectra.T)
    L, S = model.fit(max_iter=10000, iter_print=100)
    cov = np.cov(L) 
    mu = np.median(spectra, axis=0)

    # Make sure we are positive definite
    evecs,evals,T  = np.linalg.svd(cov)
    evals[evals<1e-6]=1e-6
    cov = evecs @ np.diag(evals) @ T
    cov = cov + np.eye(rows) 

    # Now we form a bad pixel map
    bad = np.zeros((rows, columns))
    for j in range(columns):
        improvements = np.zeros((len(frames),rows))
        print('forming bad pixel map:',j,'/',columns)
        for k in range(len(frames)):
            x = frames[k][j,:] 
            npca = args.npca
            xproj = (x - mu)[np.newaxis,:] @ evecs[:,:npca]
            xproj = xproj @ (evecs[:,:npca]).T + mu
            improvements[k,:] = abs(xproj - x)
        improvement = np.median(improvements,axis=0)

        # Threshold the improvements
        median_improvement = np.median(improvement)
        stdev_improvement = np.std(improvement)
        bad[:,j] = improvement > median_improvement + (stdev_improvement * args.thresh)
        print(sum(bad[:,j]),np.where(bad[:,j])[0])

    # Write the bad pixel map
    with open(args.output_bad,'wb') as fbout:
       np.array(bad, dtype=np.float32).tofile(fbout)

    # Finally, fix all the frames
    with open(args.output,'wb') as fout:
       with open(args.input,'rb') as fin:
          for line in range(lines):
              # Read a frame of data
              if line%10==0:
                  logging.info('Line '+str(line))
                  print(line)
              frame = np.fromfile(fin, count=nframe, dtype=dtype)
              frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)

              for i in range(frame.shape[1]):

                   tofix = bad[:,i]
                   nottofix = np.logical_not(tofix)
                   tofix=np.where(tofix)[0]
                   nottofix=np.where(nottofix)[0]

                   a,b = conditional_gaussian(mu, cov, tofix, nottofix, frame[nottofix,i])
                   plt.plot(frame[:,i],'r')
                   frame[tofix,i] = a
                   plt.plot(frame[:,i],'b')
                   plt.show()
              np.array(frame,dtype=np.float32).tofile(fout)
         
    print('done') 

if __name__ == '__main__':

    main()
