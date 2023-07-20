import numpy as np
from itertools import combinations
from glob import glob
from spt3g.mapspectra import curved_sky as cs
import os
import healpy as hp
from spt3g.util import healpix_tools as hpt

root = '/sptgrid/user/rgualtie/'
wmap_names = glob(root+'wmapMock/wmap_*/*GHz/coadd_wmap_*.g3')
planckLFI_names = glob(root+'PlanckMock/LFI_SkyMap_*/*GHz/coadd_LFI_SkyMap_*.g3')
planckHFI_names = glob(root+'PlanckMock/HFI_SkyMap_*/*GHz/coadd_HFI_SkyMap_*.g3')
spt3g_names = glob(root+'spt3g_full_depth_maps/*')

crosses = []
for comb in combinations(wmap_names+planckLFI_names+planckHFI_names+spt3g_names, 2):
    if os.path.isfile('/home/rgualtie/experiments_cross_spectra/'+os.path.splitext(os.path.basename(comb[0]))[0]+'_cross_'+os.path.splitext(os.path.basename(comb[1]))[0]+'.npy'):
        continue
    map1 = [fr for fr in core.G3File(comb[0])][0];maps.RemoveWeights(map1, zero_nans=True)
    map2 = [fr for fr in core.G3File(comb[1])][0];maps.RemoveWeights(map2, zero_nans=True)
    nside1 = hp.npix2nside(np.asarray(map1['T']).shape[0])
    nside2 = hp.npix2nside(np.asarray(map2['T']).shape[0])
    if nside1>nside2:
        for k in map1.keys():
            temp = map1.pop(k)
            temp = hp.ud_grade(temp, nside2)
            map1[k] = temp
            del temp
        mask = hp.ud_grade(mask, nside2)
    if nside1<nside2:
        for k in map2.keys():
            temp = map2.pop(k)
            temp = hp.ud_grade(temp, nside1)
            map2[k] = temp
            del temp
        mask = hp.ud_grade(mask, nside1)        
    
    ellb, dlsb = cs.spectrum_spice(map1, map2, lmin=8, lmax=6143, bin_width=25, mask=mask, return_error=False, lfac=True, tolerance=1.e-7)
    np.save('/home/rgualtie/experiments_cross_spectra/'+os.path.splitext(os.path.basename(comb[0]))[0]+'_cross_'+os.path.splitext(os.path.basename(comb[1]))[0]+'.npy', dlsb)
    del map1
    del map2 

