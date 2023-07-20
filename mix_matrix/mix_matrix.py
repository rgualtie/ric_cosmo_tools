#mix matrix
import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import ric_tools as rt
import subprocess
from spt3g import core, mapmaker, util, mapspectra, maps, std_processing

parser = argparse.ArgumentParser()
parser.add_argument('--lmax', type=int, default=300,
                   help = 'Maximum ell, default=300')
parser.add_argument('--nside', type=int, default=1024,
                   help = 'Nside of the output maps, default=1024')
parser.add_argument('--path', type=str, default='/big_scratch/rgualtie/mix_matrix/',
                   help = 'Path to the working directory')
parser.add_argument('--nobs', type=int, default=512,
                   help = 'Number of ObsIds for the mock-observation, default=512')
parser.add_argument('--nsims', type=int, default=100,
                   help = 'Number of sims per spectrum, default=100')
parser.add_argument('--debug',action = 'store_true',
                   help = 'Run in debug mode, overwrite nobs=64, nsims=1')
pargs = parser.parse_args()

if not os.path.exists(pargs.path):
    os.makedirs(pargs.path)

# Commodities
mask = maps.fitsio.load_skymap_fits('/sptgrid/user/rgualtie/ApodMaskPS2048.fits')
mask = np.asarray(mask['T'])
mask[mask==0]=np.nan
mask = hp.ud_grade(mask, pargs.nside)

th_specs = np.genfromtxt('/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing/base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat', unpack=True)
th_specs = np.vstack((th_specs, np.zeros_like(th_specs[0]), np.zeros_like(th_specs[0])))
lmin = 8
lmax = pargs.lmax
lbin = 25
th_specs = th_specs[1:]
th_ellb, th_spec_bin = rt.bin_spectrum(th_specs, lmin=lmin, lmax=None, binwidth=lbin, return_error=False)

#-----------------------------------------------------------------------------------------------------------------
files = glob(pargs.path+'spectra_*.npy')
# Check if the spectra are there or create them
if not files:
    lmin=lmin
    lmax=pargs.lmax
    bin_width=lbin
    start_bins = np.arange(lmin,lmax,bin_width)
    titles = ['TT','EE','BB','TE','EB','TB']
    temp = np.zeros_like(th_specs)
    ell = np.arange(th_specs.shape[1])
    for i in range(6):
        for b in start_bins:
            temp[i][b: b+bin_width] = th_specs[i][b: b+bin_width]
            with open(pargs.path+'spectra_'+titles[i]+'_'+str(b).zfill(4)+'.npy', 'wb') as f:
                np.save(f, temp)
            temp = np.zeros_like(th_specs)
#-----------------------------------------------------------------------------------------------------------------

nobs = pargs.nobs
# Number of sims for each power spectrum
nsims = pargs.nsims
# Run for one EE bin only and one sim
if pargs.debug:
    print('Debug mode ON')
    nobs = 512
    nsims = 1
    files = [pargs.path+'spectra_EE_0033.npy']
    
# Folder on the grid is needed to mock-observe
if not os.path.exists('/sptgrid/user/rgualtie/mix_matrix'):
    os.system('gfal-mkdir -p gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/')

nfiles = len(files)
#make the maps
for f in files:
    cls = np.load(f)
    for i in range(nsims):
        sd = 1000+i
        np.random.seed(sd)
        mm = hp.synfast(cls, nside=pargs.nside, pol=True, new=True)
        mapname = f.split('/')[-1][:-4]+'_s'+str(sd).zfill(4)#+'.fits'
        hp.write_map(pargs.path+mapname+'.fits', m=mm, coord='G', overwrite=True)
        os.system('gfal-copy '+pargs.path+mapname+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix')
        os.system('rm '+pargs.path+mapname+'.fits')
        #mock-observe each map
        os.system('for f in 90 150 220;do python /home/rgualtie/submit_mock_obs_map_dag.py --submit --inmap /sptgrid/user/rgualtie/mix_matrix/'+mapname+'.fits --healpix --nside '+str(pargs.nside)+' --band $f --randomobs --nobs '+str(nobs)+' --mixmat;done')
        # Wait until the mock is done
        print('Mock-observing %s.fits, hold on..'%mapname)
        time.sleep(1)
        # How do I check that all is ready?
        # Here the code should fork and open threads so that while the mock is in progress for one map the others are submitted
        # Need some improvement
        while 1:
            os.system('condor_q -nobatch > '+pargs.path+'temp.txt')
            queue = np.genfromtxt(pargs.path+'temp.txt', unpack=True, dtype=str, delimiter='\t') 
            queue = queue[2:-3]
            running_procs = []
            for i in range(len(queue)):
                running_procs.append(queue[i].split()[-1])
            if any(mapname in s for s in running_procs):
                # wait 10min and check again
                print('The grid is still processing the map..check in 10min')
                os.system('rm '+pargs.path+'temp.txt')
                time.sleep(600)
            else:
                print('%s is ready for coadd'%mapname)
                break
        # Coadd the mock observed maps, generate the spectra and cleanup
        if pargs.debug:
            os.system('python /home/rgualtie/src/coadder.py --inmap '+mapname+' --chunk_len 64 --spectra --mixmat --cleanup --debug')
        else:
            os.system('python /home/rgualtie/src/coadder.py --inmap '+mapname+' --chunk_len 64 --spectra --mixmat --cleanup')
            os.system('gfal-rm gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/'+pargs.path+mapname+'.fits')
