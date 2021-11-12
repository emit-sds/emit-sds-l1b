from scipy.signal import find_peaks 
from scipy.stats import norm
from astropy import modeling
from scipy.interpolate import interp1d
import json

X = envi.open('frames.hdr').load()

lines, samps, bands = X.shape
half = int(samps/2)
ghost = np.zeros((480,480))
ghost_spatial = np.zeros((1280,1280))
center = 646.5

def find_peaks(x, thresh=2, contrast_thresh = 0.9):
    
    # Find channels exceeding one standard deviation of the mean
    n = len(x)
    stdev = x.std()
    divergence = x-x.mean()
    indices = np.where(divergence>(thresh*stdev))[0]
    
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
        magnitudes = divergence[group_indices]
        
        if len(group_indices)<3:
            total = magnitudes.sum()
            avg_index = np.sum(group_indices * magnitudes / total)
            peaks.append((avg_index,max(magnitudes),np.std(group_indices)))
            continue
        fitter = modeling.fitting.LevMarLSQFitter()

        model = modeling.models.Gaussian1D(amplitude=np.max(magnitudes),
                                           mean=np.mean(group_indices),
                                           stddev=np.std(group_indices))   # depending on the data you need to give some initial values
        fitted_model = fitter(model, group_indices, magnitudes)
        
        if fitted_model.mean[0] > 1279 or fitted_model.mean[0] < 0:
            continue
        
        peaks.append((fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]))
    max_amplitude = max([peaks[j][1] for j in range(len(peaks))])
    
    for i in range(len(peaks)-1,-1,-1):
        if peaks[i][1] < max_amplitude:
            del peaks[i]
    return np.array(peaks)[0]

pairs = []
for i in range(lines):
    
    x = np.squeeze(X[i,:,:]).T
    x = gaussian_filter(x,[1,1])
    L,R = x[:,:half], x[:,half:]
    
    rt_max_colwise  = np.max(R, axis=0)
    lft_max_colwise = np.max(L,axis=0)
    rt_max_rowwise  = np.max(R, axis=1)
    lft_max_rowwise = np.max(L,axis=1)
    rt_sum = np.sum(R[R > 0.0001])
    lft_sum = np.sum(L[L> 0.0001])
   
    best_left_spectral  = find_peaks(lft_max_rowwise)
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
blurs_spatial = interp1d([x[0] for x in pairs], [x[3] for x in pairs],
                         fill_value=np.mean([x[3] for x in pairs]),
                         bounds_error=False)(range(480))
blurs_spectral = interp1d([x[0] for x in pairs], [x[4] for x in pairs],
                         fill_value=np.mean([x[4] for x in pairs]),
                         bounds_error=False)(range(480))

lines, extents, blurs_spatial, blurs_spectral, intensities = [],[],[],[],[]
dataset = []
for i in range(len(pairs)):
    print(pairs[i])
    # recognize a gap in the lines
    if len(dataset)>0 and pairs[i][1]>dataset[-1][1] or i==(len(pairs)-1):
        
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
        blurs_spatial.append(np.mean([d[3] for d in dataset]))
        blurs_spectral.append(np.mean([d[4] for d in dataset]))
        dataset = [pairs[i]]
    else:
        dataset.append(pairs[i])
        
blur_spatial = np.zeros(480)
blur_spectral = np.zeros(480)
intens = np.zeros(480)
ghost = np.zeros((480,480))
for i in range(480):
    for j,extent in enumerate(extents):
        if i>=extent[0] and i<extent[1]:
            target = int(round(lines[j][0] * i + lines[j][1]))
            ghost[i,target] = 1
            blur_spatial[i] = blurs_spatial[j]
            blur_spectral[i] = blurs_spectral[j]
            intens[i] = intensities[j]
            
ghost_config = {'center':center, 'orders':[]}
for j,extent in enumerate(extents):
    ghost_config['orders'].append({'blur_spatial':float(blurs_spatial[j]),
                                  'blur_spectral':float(blurs_spectral[j]),
                                  'slope':float(lines[j][0]),
                                  'offset':float(lines[j][1]),
                                  'extent':(int(extents[j][0]),int(extents[j][1])),
                                  'intensity':float(intensities[j])})

with open('emit_ghost.json','w') as fout:
    json.dump(ghost_config,fout)

