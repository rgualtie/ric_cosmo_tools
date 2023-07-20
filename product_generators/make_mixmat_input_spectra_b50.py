import numpy as np
import os, time, argparse
from glob import glob
import ric_tools as rt

parser = argparse.ArgumentParser()
parser.add_argument('--freq', type=str, default='90',
                   help = 'Frequency, adds an MBB dust spectrum to the theoretical spectrum')
parser.add_argument('--path', type=str, default='/big_scratch/rgualtie/mix_matrix/input_spectra_b50/',
                   help = 'Path to the working directory')
pargs = parser.parse_args()

if not os.path.exists(pargs.path):
    os.makedirs(pargs.path)

th_specs = np.genfromtxt('/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat', unpack=True)
th_specs = np.vstack((th_specs, np.zeros_like(th_specs[0]), np.zeros_like(th_specs[0])))
lmin = 8
lbin = 50
th_specs = th_specs[1:]
dust = np.load('spectra/'+pargs.freq+'GHz/'+pargs.freq+'GHz_MBB_spectrum.npy')
inspec = np.zeros_like(th_specs[:,:len(dust)])
for i in range(6):
    inspec[i] = th_specs[i, :len(dust)]+dust
#inspec_bin = rt.bin_spectrum(inspec, lmin=lmin, lmax=None, binwidth=lbin, return_error=False)[1]

#----------------------------------------------------------------------------------------------------------------
files = glob(pargs.path+'spectra_*.npy')
# Check if the spectra are there or create them
if not files:
    lmin=lmin
    lmax=6143
    bin_width=lbin
    start_bins = np.arange(lmin, lmax+bin_width, bin_width)
    titles = ['TT','EE','BB','TE','EB','TB']
    temp = np.zeros_like(inspec)
    ell = np.arange(inspec.shape[1])
    for i in range(6):
        for b in start_bins:
            temp[i][b: b+bin_width] = inspec[i][b: b+bin_width]
            with open(pargs.path+'spectra_'+titles[i]+'_'+str(b).zfill(4)+'.npy', 'wb') as f:
                np.save(f, temp)
            temp = np.zeros_like(inspec)
#----------------------------------------------------------------------------------------------------------------

# Move the spectra to the grid
os.system('gfal-copy -r /big_scratch/rgualtie/mix_matrix/input_spectra_b50/ gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/input_spectra_b50/'+pargs.freq+'/')

