import numpy as np
from glob import glob

def avdict(dictionary):
    res = 0
    for val in dictionary.values():
        res += val
    return res/len(dictionary)

def stddict(dictionary):
    listr = []
    for val in dictionary.values():
        listr.append(val)
    return np.std(listr, axis=0)

freq = '90'
spectra = {'TT':'T','EE':'E','BB':'B','TE':'TE','EB':'EB','TB':'TB'}
for i,s in enumerate(spectra.keys()):
    specs = []
    for j in range(500):
        try:
            specs.append(np.load('/sptgrid/user/rgualtie/EnoB/output_spectra_bkmask/'+freq+'GHz/spectra_EnoBlens_s1'+str(j).zfill(3)+'_idx'+str(j).zfill(3)+'_'+freq+'GHz_bkmask.npz')['dlsb'][i])
        except:
            continue
    specs = np.array(specs)
    cov = np.cov(specs.T)

    noise_f = np.sort(glob('/sptgrid/user/rgualtie/spectra/'+freq+'GHz/noise_bkmask/'+freq+'GHz_noise_spectra_bkmask_0*.npz'))
    nspectra = {}
    nspectra['T']={};nspectra['E']={};nspectra['B']={}
    nspectra['TE']={};nspectra['EB']={};nspectra['TB']={}
    for k,f in enumerate(noise_f):
        if 'Cl' in f:
            continue
        nspectra['T'][k] =  np.load(f)['dlsb'][0]
        nspectra['E'][k] = np.load(f)['dlsb'][1]
        nspectra['B'][k] = np.load(f)['dlsb'][2]
        nspectra['TE'][k] = np.load(f)['dlsb'][3]
        nspectra['EB'][k] = np.load(f)['dlsb'][4]
        nspectra['TB'][k] = np.load(f)['dlsb'][5]

    nTTmean = avdict(nspectra['T'])
    nEEmean = avdict(nspectra['E'])
    nBBmean = avdict(nspectra['B'])
    nTEmean = avdict(nspectra['TE'])
    nEBmean = avdict(nspectra['EB'])
    nTBmean = avdict(nspectra['TB'])

    nTTstd = stddict(nspectra['T'])
    nEEstd = stddict(nspectra['E'])
    nBBstd = stddict(nspectra['B'])
    nTEstd = stddict(nspectra['TE'])
    nEBstd = stddict(nspectra['EB'])
    nTBstd = stddict(nspectra['TB'])

    noise_spectra = np.array([nTTmean, nEEmean, nBBmean, nTEmean, nEBmean, nTBmean])
    noise_std = np.array([nTTstd, nEEstd, nBBstd, nTEstd, nEBstd, nTBstd])
    nspecs = []
    for v in nspectra[spectra[s]].keys():
        nspecs.append(1e6*nspectra[spectra[s]][v])
    nspecs = np.array(nspecs)
    ncov = np.cov(nspecs.T)
    np.save(freq+'GHz_'+s+'_SpNCov_bkmask.npy',cov+ncov)
