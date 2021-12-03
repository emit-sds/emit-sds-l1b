from scipy.signal import find_peaks 
from scipy.stats import norm
from astropy import modeling
import sys
import numpy as np
from spectral.io import envi
from scipy.interpolate import interp1d
import json
from scipy.ndimage import gaussian_filter
np.set_printoptions(precision=5)


X = envi.open(sys.argv[1]+'.hdr').load()

lines, samps, bands = X.shape
half = int(samps/2)
center = 649

def find_peaks(xo, thresh = 20.0, min_DN = 5):
    
    x = gaussian_filter(xo[np.newaxis,:],[1,1])[0,:]

    # Find channels exceeding four standard deviation of the mean
    n = len(x)
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
        
        if fitted_model.amplitude[0] < min_DN:
            continue
           
        peaks.append((fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]))

    return np.array(peaks)

pairs = []

for peak in [0,1]:
  
  for i in range(lines):
    
    x = np.squeeze(X[i,:,:]).T
    ghosts = find_peaks(x[:,359])
    sources = find_peaks(x[:,939])
 
    if len(ghosts)>peak:

        ghost = ghosts[peak]
        source = np.array(sources[0])
        print(source[0],ghost[0])
        source_loc = source[0]
        ghost_loc = ghost[0]
        source_energy = source[1]
        ghost_energy = ghost[1]
        ghost_stdev = ghost[2]
        
        pairs.append((source_loc,
                      ghost_loc,
                      ghost_energy/source_energy,
                      ghost_stdev))

#sys.exit(0)

if False:
    
    rt_max_colwise  = np.max(R, axis=0)
    lft_max_colwise = np.max(L,axis=0)
    rt_max_rowwise  = np.max(R, axis=1)
    lft_max_rowwise = np.max(L,axis=1)
    rt_sum = np.sum(R[R > 0.0001])
    lft_sum = np.sum(L[L> 0.0001])
   
    best_right_spectral = find_peaks(rt_max_rowwise)
    best_left_spatial   = find_peaks(lft_max_colwise)
    best_right_spatial  = find_peaks(rt_max_colwise)
    
    if L.max()>R.max():
        source_loc_spatial, source_amplitude_spatial, source_stdev_spatial = best_left_spatial
        target_loc_spatial, target_amplitude_spatial, target_stdev_spatial = best_right_spatial
        source_loc_spectral, source_amplitude_spectral, source_stdev_spectral = best_left_spectral
        target_loc_spectral, target_amplitude_spectral, target_stdev_spectral = best_right_spectral
        source_energy, target_energy = lft_sum, rt_sum
    else:
        source_loc_spatial, source_amplitude_spatial, source_stdev_spatial = best_right_spatial
        target_loc_spatial, target_amplitude_spatial, target_stdev_spatial = best_left_spatial
        source_loc_spectral, source_amplitude_spectral, source_stdev_spectral = best_right_spectral
        target_loc_spectral, target_amplitude_spectral, target_stdev_spectral = best_left_spectral
        source_energy, target_energy = rt_sum, lft_sum
    if target_stdev_spectral > 10 or target_stdev_spatial > 10 or \
        target_stdev_spectral < 1 or target_stdev_spatial < 1:
        print('skipping target ',target_loc_spatial)
        continue
    
    pairs.append((source_loc_spectral,
                  target_loc_spectral,
                  target_energy/source_energy,
                  target_stdev_spatial,
                  target_stdev_spectral))

pairs.sort()
blurs = interp1d([x[0] for x in pairs], [x[3] for x in pairs],
                         fill_value=np.mean([x[3] for x in pairs]),
                         bounds_error=False)(range(480))

lines, extents, blurs, intensities = [],[],[],[]
dataset = []
for i in range(len(pairs)):
    print(pairs[i])

    source_loc, ghost_loc, intens, stdev = pairs[i]

    # recognize a gap in the data sequence, create a new line
    if len(dataset)>0 and ghost_loc>dataset[-1][1] or i==(len(pairs)-1):
        
        # ignore really small segments
        if len(dataset)<2:
            dataset = [pairs[i]]
            continue
            
        # fit a line
        slope, offset = np.polyfit([x[0] for x in dataset],
                                   [y[1] for y in dataset],1)
        lines.append((slope, offset)) 
        extents.append((int(round(dataset[0][0])),
                        int(round(pairs[i][0]))))
        intensities.append(np.mean([d[2] for d in dataset]))
        blurs.append(np.mean([d[3] for d in dataset]))
        dataset = [pairs[i]]

    else:
        dataset.append(pairs[i])
        
if False:
    blur = np.zeros(480)
    intens = np.zeros(480)
    ghost = np.zeros((480,480))
    for i in range(480):
        for j,extent in enumerate(extents):
            if i>=extent[0] and i<extent[1]:
                target = int(round(lines[j][0] * i + lines[j][1]))
                ghost[i,target] = 1
                blur[i] = blurs[j]
                intens[i] = intensities[j]
            
ghost_config = {'center':center, 'orders':[]}
for j,extent in enumerate(extents):
    ghost_config['orders'].append({'blurs':float(blurs[j]),
                                  'slope':float(lines[j][0]),
                                  'offset':float(lines[j][1]),
                                  'extent':(int(extents[j][0]),int(extents[j][1])),
                                  'intensity':float(intensities[j])})

with open(sys.argv[2],'w') as fout:
    s=json.dumps(ghost_config,indent=2)
    fout.write(s)

