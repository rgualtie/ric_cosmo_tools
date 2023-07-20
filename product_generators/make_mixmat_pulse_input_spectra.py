import matplotlib.pylab as plt
import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import ric_tools as rt
import os
from spt3g import core, mapmaker, util, mapspectra, maps, std_processing

parser = argparse.ArgumentParser("")
parser.add_argument('--lmax', type=int, default=500,
                   help = 'Maximum ell, default=500')
parser.add_argument('--bin_step', type=int, default=25)
parser.add_argument('--path', type=str, default='/big_scratch/rgualtie/mix_matrix/input_pulse_spectra/',
                   help = 'Path to the working directory')
pargs = parser.parse_args()

if not os.path.exists(pargs.path):
    os.makedirs(pargs.path)

th_specs = np.genfromtxt('/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat', unpack=True)
th_specs = np.vstack((th_specs, np.zeros_like(th_specs[0]), np.zeros_like(th_specs[0])))
lmin = 8
lmax = pargs.lmax
lbin = 1
th_specs = th_specs[1:]
th_ellb, th_spec_bin = rt.bin_spectrum(th_specs, lmin=lmin, lmax=None, binwidth=lbin, return_error=False)

#----------------------------------------------------------------------------------------------------------------
store = True
plot = True

if store:
    #files = glob(pargs.path+'pulse_spectra_*.npy')
    # Check if the spectra are there or create them
    #if not files:
    lmin=lmin
    lmax=pargs.lmax
    bin_width=pargs.bin_step
    start_bins = np.arange(lmin, lmax+bin_width, bin_width)
    titles = ['TT','EE','BB','TE','EB','TB']
    temp = np.zeros_like(th_specs)
    ell = np.arange(th_specs.shape[1])
    fig,ax = plt.subplots(1,3,sharex=True)
    plt.suptitle('Input pulses spectra')
    for i in range(3):
        for b in start_bins:
            temp[i][b: b+lbin] = th_specs[i][b: b+lbin]
            with open(pargs.path+'pulse_spectra_'+titles[i]+'_'+str(b).zfill(4)+'.npy', 'wb') as f:
                np.save(f, temp)
            if plot:
                ax[i].plot(ell, temp[i])
                if b==lmin:
                    ax[i].plot(th_specs[i],'k-', label='$\Lambda$-CDM')
                ax[i].set_xlim(0,500)
                ax[i].set_title(titles[i])
                ax[i].set_xlabel('Multipole $\ell$')
                if i==0:
                    ax[i].set_ylabel('$D_\ell\quad [\mu K^2]$')
                ax[i].grid();ax[i].legend()
            temp = np.zeros_like(th_specs)
#----------------------------------------------------------------------------------------------------------------

# Move the spectra to the grid
files = glob(pargs.path+'pulse_spectra_*.npy')
if not files:
    os.system('gfal-copy -r /big_scratch/rgualtie/mix_matrix/input_pulse_spectra/ gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/input_pulse_spectra/')

