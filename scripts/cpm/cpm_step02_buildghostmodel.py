# David R Thompson
import json
import numpy as np
import pylab as plt


ghost_config = {'center':322.5, 'scale_all':True, 'orders':[],
          'psf_zones':[{'extent':[20,268], 'psfs': [{'sigma':1, 'peak':0.000001},
                                                   {'sigma':50,'peak':0.01}]},
                       {'extent':[268,394], 'psfs': [{'sigma':1, 'peak':0.000001},
                                                   {'sigma':50,'peak':0.01}]}]}
          

plot = False
write = True 

with open('../../data/cpm/ghost_pointwise_edit.txt','r') as fin:
   sources, targets, intensities = [],[],[]
   for line in fin.readlines():
     if line[0] == '#':
         continue
     elif len(line)>3:
         toks = line.split()
         sources.append(float(toks[0]))
         targets.append(float(toks[1]))
         intensities.append(float(toks[-1]))
     else:
         slope, offset = np.polyfit(sources,targets,1)
         islope, ioffset = np.polyfit(sources,intensities,1)
         order = {'slope':slope, 'offset':offset, 
            'extent':(min(sources),max(sources)),
            'scaling':1.0,
            'intensity_slope':islope,
            'intensity_offset':ioffset}
         ghost_config['orders'].append(order)
         x = np.arange(min(sources),max(sources)+1)
         if plot:
             plt.figure(0)
             plt.plot(x,x*slope+offset,'r')
             plt.plot(sources,targets,'k.')
             plt.figure()
             plt.plot(sources,intensities)
         sources, targets, intensities = [],[],[]
if plot:
    plt.show()
    plt.savefig('ghost_image.png')

if write:
    with open('../../data/cpm/cpm_ghost.json','w') as fout:
        fout.write(json.dumps(ghost_config, indent=4))
         
   
  
