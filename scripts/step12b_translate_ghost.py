# David R Thompson
import numpy as np
import pylab as plt
import os, sys, glob
import pylab as plt

x = np.loadtxt('../data/ghost_pointwise.txt')
y = []
for i in range(x.shape[0]):
    n = x[i].copy()
    if len(y) == 0:
        y.append(n)
        continue
    if n[1]<y[-1][1]:
        y.append(n)
        continue
    else:
       gap = max(int(n[1])-int(y[-1][1]),1)
       for j in range(int(y[-1][1])+1,int(n[1])):
         m = n.copy()
         m[3] = m[3] / float(gap)
         m[1] = j
         y.append(m)
       n[3] = n[3] / float(gap)
       y.append(n)
       plt.plot(n[0],n[1],'ko')
plt.show()
with open('../data/ghost_pointwise_interp.txt','w') as fout:
   for n in y:
       fout.write('%i %i %10.8f %10.8f\n' % (n[0],n[1],n[2],n[3]))

