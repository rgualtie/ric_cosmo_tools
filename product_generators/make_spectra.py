import numpy as np
import os, sys, time
import argparse
import healpy as hp
from glob import glob
from spt3g import core, std_processing, maps
from spt3g.mapspectra import curved_sky as cs

parser = argparse.ArgumentParser()
parser.add_argument('--freq',type=str, default='90',
                   help = 'Frequency, default = 90')
parser.add_argument('--freq2',type=str, default='150',
                   help = 'Frequency, default = 150')
parser.add_argument('--lmin', type=int, default=8,
                   help = 'Lowest multipole, default = 8')
parser.add_argument('--lmax', type=int, default=6143,
                   help = 'Highest multipole, default = 6143')
parser.add_argument('--binw', type=int, default=25,
                   help = 'Bin width, default = 25')
parser.add_argument('--inmap',type=str,
                    help = 'input maps path on the grid, full path is not required')
parser.add_argument('--inmap2',type=str, 
                    help = 'second input maps path on the grid, full path is not required')
parser.add_argument('--mask',type=str, default='/sptgrid/user/rgualtie/ApodMaskPS2048.fits',
                    help = 'mask path on the grid, full path is required')
parser.add_argument('--output', type=str, default='spectra.npz',
                    help = 'Output name for the file, if not specified it defaults on the input map(s) name(s)')
pargs = parser.parse_args()

f_dict = {'90':0,'150':1,'220':2}

d = [fr for fr in core.G3File(pargs.inmap)][f_dict[pargs.freq]]
maps.RemoveWeights(d, zero_nans=True)
#umap = d.pop('U')
#umap *= -1
#d['U'] = umap

if pargs.inmap2:
    d2 = [fr for fr in core.G3File(pargs.inmap2)][f_dict[pargs.freq2]]
    maps.RemoveWeights(d2, zero_nans=True)
    #umap2 = d2.pop('U')
    #umap2 *= -1
    #d2['U'] = umap2
else:
    d2 = None

mask = maps.fitsio.load_skymap_fits(pargs.mask)
try:
    mask = np.asarray(mask['T'])
except:
    mask = np.asarray(mask['AP_MASK'])
mask[mask==0]=np.nan
mask = hp.ud_grade(mask, hp.npix2nside(len(np.asarray(d['T']))))

ellb, dlsb, errb = cs.spectrum_spice(map1=d, map2=d2, lmin=pargs.lmin, lmax=pargs.lmax, bin_width=pargs.binw, mask=mask, lfac=True, return_error=True, tolerance=1.e-7)
#store the spectra
if not pargs.output:
    if pargs.inmap2:
        pargs.output = pargs.inmap+'_'+pargs.freq+'GHz_cross_'+pargs.inmap2+'_'+pargs.freq2+'GHz.npz'
    else:
        pargs.output = pargs.inmap+'_'+pargs.freq+'GHz.npz'

np.savez(pargs.output, ellb=ellb, dlsb=dlsb, errb=errb)
