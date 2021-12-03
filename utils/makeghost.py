from scipy.signal import find_peaks 
from scipy.stats import norm
from astropy import modeling
import sys
import numpy as np
from spectral.io import envi
from scipy.interpolate import interp1d
import json
from scipy.ndimage import gaussian_filter
import pylab as plt
from skimage.measure import LineModelND, ransac
np.set_printoptions(precision=5)


X = envi.open(sys.argv[1]+'.hdr').load()

lines, samps, bands = X.shape
half = int(samps/2)
center = 649

def find_peaks(x, xo, thresh = 20.0):
    
    # Find channels exceeding four standard deviation of the mean
    n = len(x)

    x[:25] = 0
    x[315:] = 0

    indices = np.where(x>thresh)[0]
    if len(indices)<1:
        return []
    
    # Separate into contiguous sequences
    groups = [[indices[0]]]
    for i in range(1,len(indices)):
        if indices[i]==(groups[-1][-1] + 1):
            groups[-1].append(indices[i])
        else:
            groups.append([indices[i]])
            
    # Fit a gaussian shape to each sequence
    peaks = []
    for group_indices in groups:
   
        # use unfiltered values
        magnitudes = xo[group_indices]
        
        if len(group_indices)<3:
            continue
        fitter = modeling.fitting.LevMarLSQFitter()

        model = modeling.models.Gaussian1D(amplitude=np.max(magnitudes),
                                           mean=np.mean(group_indices),
                                           stddev=np.std(group_indices))   # depending on the data you need to give some initial values
        fitted_model = fitter(model, group_indices, magnitudes)
        
        if fitted_model.mean[0] > 1279 or fitted_model.mean[0] < 0:
            continue
           
        peaks.append((fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]))

    return np.array(peaks)

pairs = []

# enforce negative slope
def model_valid(model, data):
  origin, direction = model.params
  slope = direction[1]/direction[0]
  if slope>0:
      return False
  return True



for peak in [0,1]:
  
  for i in range(lines):
    
    x = np.squeeze(X[i,:,:]).T
    xf = gaussian_filter(x,[2,2])
    ghosts = find_peaks(x[:,359],xf[:,359])
    sources = find_peaks(x[:,939],xf[:,939])
 
    if len(ghosts)>peak:

        ghost = ghosts[peak]
        source = np.array(sources[0])
        source_loc = source[0]
        ghost_loc = ghost[0]
        source_energy = source[1]
        ghost_energy = ghost[1]
        ghost_stdev = ghost[2]

        # filter out a bad data segment
        if source_loc > 258 and source_loc < 264:
            continue
        
        pairs.append((source_loc,
                      ghost_loc,
                      ghost_energy/source_energy,
                      ghost_stdev))

#pairs.sort()
#blurs = interp1d([x[0] for x in pairs], [x[3] for x in pairs],
#                        fill_value=np.mean([x[3] for x in pairs]),
#                        bounds_error=False)(range(480))



ghost_config = {'center':center, 'orders':[]}
data = np.array([[s[0],s[1]] for s in pairs])
intensities = np.array([s[2] for s in pairs])
blurs =  np.array([s[3] for s in pairs])

for i in range(10):
  # fit line using all data
  model = LineModelND()
  model.estimate(data)

  # robustly fit line only using inlier data with RANSAC algorithm
  try:
      model_robust, inliers = ransac(data, LineModelND, 
                                     is_model_valid = model_valid,
                                     min_samples=4,
                                     residual_threshold=2, max_trials=10000)
  except ValueError:
      break
  if model_robust is None:
      break

  # Extract the parameters
  origin, direction = model_robust.params
  slope = direction[1]/direction[0]
  offset = origin[1]-slope*origin[0]
  extent = (min(data[inliers,0]),max(data[inliers,0]))
  blur = np.mean(blurs[inliers])
  intensity = np.mean(intensities[inliers])

  # Remove the datapoints from the line we have fit
  new_indices = np.ones(len(data))>0
  new_indices[inliers] = False
  data = data[new_indices,:]
  intensities = intensities[new_indices]
  blurs = blurs[new_indices]
  
  x = np.arange(extent[0],extent[1]+1)
  plt.plot(x,x*slope+offset,'r-')

  ghost_config['orders'].append({'blur':float(blur),
                                  'slope':float(slope),
                                  'offset':float(offset),
                                  'extent':(int(extent[0]),int(extent[1])),
                                  'intensity':float(intensity)})

with open(sys.argv[2],'w') as fout:
    s=json.dumps(ghost_config,indent=2)
    fout.write(s)
          
plt.plot([x[0] for x in pairs],[x[1] for x in pairs],'k.')
plt.show()


