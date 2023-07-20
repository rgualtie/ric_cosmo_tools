import numpy as np
import os, sys, time
import argparse
import healpy as hp
from glob import glob
from spt3g import core, std_processing, maps
from spt3g.mapspectra import curved_sky as cs

parser = argparse.ArgumentParser()
parser.add_argument('--lmin', type=int, default=8,
                   help = 'Lowest multipole, default = 8')
parser.add_argument('--lmax', type=int, default=6143,
                   help = 'Highest multipole, default = 6143')
parser.add_argument('--binw', type=int, default=25,
                   help = 'Bin width, default = 25')
parser.add_argument('--inmap',type=str,
                    help = 'input maps path on the grid, full path is not required')
#parser.add_argument('--inmap2',type=str, 
#                    help = 'second input maps path on the grid, full path is not required')
parser.add_argument('--mask',type=str, default='/sptgrid/user/rgualtie/ApodMaskPS2048.fits',
                    help = 'mask path on the grid, full path is required')
parser.add_argument('--output', type=str, default='spectra.npz',
                    help = 'Output name for the file, if not specified it defaults on the input map(s) name(s)')
parser.add_argument('--lfac', action = 'store_true',
                   help = 'If true calculate Dls instead of Cls')
parser.add_argument('--auto', action = 'store_true')
parser.add_argument('--cross', action = 'store_true')
parser.add_argument('--freq', type=str, default='90')
parser.add_argument('--freq2', type=str, default='150')

pargs = parser.parse_args()

freqs = {'90':0,'150':1,'220':2}

d = [fr for fr in core.G3File(pargs.inmap)]
for i in range(3):
    maps.RemoveWeights(d[0], zero_nans=True)

mask = maps.fitsio.load_skymap_fits(pargs.mask)
try:
    mask = np.asarray(mask['T'])
except:
    mask = np.asarray(mask['AP_MASK'])
mask[mask==0]=np.nan
mask = hp.ud_grade(mask, hp.npix2nside(len(np.asarray(d[0]['T']))))
if pargs.auto:
    ellb, dlsb, errb = cs.spectrum_spice(map1=d[freqs[pargs.freq]], map2=None, lmin=pargs.lmin, lmax=pargs.lmax, bin_width=pargs.binw, mask=mask, lfac=pargs.lfac, return_error=True, tolerance=1.e-7)
if pargs.cross:
    ellb, dlsb, errb = cs.spectrum_spice(map1=d[freqs[pargs.freq]], map2=d[freqs[pargs.freq2]], lmin=pargs.lmin, lmax=pargs.lmax, bin_width=pargs.binw, mask=mask, lfac=pargs.lfac, return_error=True, tolerance=1.e-7)

np.savez(pargs.output, ellb=ellb, dlsb=dlsb, errb=errb)
