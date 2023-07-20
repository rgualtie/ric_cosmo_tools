import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import ric_tools as rt
import os
from spt3g import core, mapmaker, util, mapspectra, maps, std_processing

parser = argparse.ArgumentParser()
parser.add_argument('--lmax', type=int, default=550,
                   help = 'Maximum ell, default=550')
parser.add_argument('--band', type=str, default=90,
                   help = 'Frequency, default=90')
parser.add_argument('--path', type=str, default='/big_scratch/rgualtie/mix_matrix/input_spectra/',
                   help = 'Path to the working directory')
parser.add_argument('--lbin', type=int, default=25,
                   help = 'delta ell bin width, default=25')
parser.add_argument('--fg', action = 'store_true',
                   help = 'If true, adds the foregrounds to the input spectra')
pargs = parser.parse_args()

if not os.path.exists(pargs.path):
    os.makedirs(pargs.path)

th_specs = np.genfromtxt('/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat', unpack=True)
th_specs = np.vstack((th_specs, np.zeros_like(th_specs[0]), np.zeros_like(th_specs[0])))
lmin = 8
lmax = pargs.lmax
lbin = pargs.lbin
th_specs = th_specs[1:, 8:6144]
if pargs.fg:
    th_specs += np.array([np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy')])
    th_specs += np.array([np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy')])
#th_ellb, th_spec_bin = rt.bin_spectrum(th_specs, lmin=lmin, lmax=None, binwidth=lbin, return_error=False)

#----------------------------------------------------------------------------------------------------------------
if pargs.fg:
    files = glob(pargs.path+pargs.band+'GHz/'+pargs.band+'GHz_spectra_*.npy')
else:
    files = glob(pargs.path+'spectra_*.npy')
# Check if the spectra are there or create them
if not files:
    lmin=lmin
    lmax=pargs.lmax
    bin_width=lbin
    start_bins = np.arange(lmin, lmax+bin_width, bin_width)
    titles = ['TT','EE','BB','TE','EB','TB']
    temp = np.zeros_like(th_specs)
    ell = np.arange(th_specs.shape[1])
    for i in range(6):
        for b in start_bins:
            temp[i][b: b+bin_width] = th_specs[i][b: b+bin_width]
            if pargs.fg:
                fname = pargs.path+pargs.band+'GHz/'+pargs.band+'GHz_spectra_'+titles[i]+'_'+str(b).zfill(4)+'_fg.npy'
            else:
                fname = pargs.path+'spectra_'+titles[i]+'_'+str(b).zfill(4)+'.npy'
            with open(fname, 'wb') as f:
                np.save(f, temp)
            temp = np.zeros_like(th_specs)
#----------------------------------------------------------------------------------------------------------------

# Move the spectra to the grid
if pargs.fg:
    os.system('gfal-copy -r /big_scratch/rgualtie/mix_matrix/input_spectra/'+pargs.band+'GHz/ gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/input_spectra/'+pargs.band+'GHz/')
else:
    os.system('gfal-copy -r /big_scratch/rgualtie/mix_matrix/input_spectra/ gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/input_spectra/')

