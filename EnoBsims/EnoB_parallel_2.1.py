# Generate EnoB+lensing simulated spectra
import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import sys
parser = argparse.ArgumentParser()
parser.add_argument('--nside', type=int, default=2048,
                   help = 'Nside of the output maps, default=2048')
parser.add_argument('--infile', type=str, default=os.environ['SPT3G_SOFTWARE_PATH']+'/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat',
                   help = 'Filename of the .npy file containing the input spectra')
parser.add_argument('--band', type=str, default='90',
                   help = 'Frequency band the mock-observation, default=90')
parser.add_argument('--group', type=str, default='EnoB',
                   help = '/sptgrid/user/$USER/<group> Subfolder where the results will be stored')
parser.add_argument('--simidx', type=int, default=0,
                   help = 'Sim index, default=0')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of files used for the mock-observation, default=256')
parser.add_argument('--oversampling', type=int, default=1,
                   help = 'Input map Nside multiplier, default=1')
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Delete all the intermediate products')
parser.add_argument('--bkmask', action = 'store_true',
                   help = 'Use bk14 mask')
parser.add_argument('--binw', type=str, default='25')
parser.add_argument('--fg', action = 'store_true',
                   help = 'If true, adds the foregrounds to the input spectra')
pargs = parser.parse_args()

# example call: $python src/EnoB_parallel_2.1.py --submit --cleanup --band 90 --simidx 0
if pargs.bkmask:
    if os.path.exists('/sptgrid/user/'+os.environ.get('USER')+'/'+pargs.group+'/output_spectra_bkmask/'+pargs.band+'GHz/spectra_'+pargs.group+'lens_s1'+str(pargs.simidx).zfill(3)+'_idx'+str(pargs.simidx).zfill(3)+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz'):
        print('output_spectrum for %i index already created, moving on' %pargs.simidx)
        sys.exit()
else:
    if os.path.exists('/sptgrid/user/'+os.environ.get('USER')+'/'+pargs.group+'/output_spectra/'+pargs.band+'GHz/spectra_'+pargs.group+'lens_s1'+str(pargs.simidx).zfill(3)+'_idx'+str(pargs.simidx).zfill(3)+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz'):
        print('output_spectrum for %i index already created, moving on' %pargs.simidx)
        sys.exit() 

f = pargs.infile
# make the maps
dls = np.genfromtxt(pargs.infile, unpack=True)[1:, 8:6144]
if pargs.fg:
    dls += np.array([np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_MBB_spectrum.npy')])
    dls += np.array([np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy'),
                     np.load(os.environ.get('HOME')+'/spectra/'+pargs.band+'GHz/'+pargs.band+'GHz_sync_spectrum.npy')])
ell = np.arange(2,len(dls[0])+2)
lfac = ell*(ell+1)/2./np.pi
cls = dls/lfac

sd = 1000+pargs.simidx
np.random.seed(sd)
mm = hp.synfast(cls, nside=pargs.oversampling*pargs.nside, lmax=3*pargs.nside-1, pol=True, new=True)

mapname = ''+pargs.group+'lens_s'+str(sd).zfill(4)+'_idx'+str(pargs.simidx).zfill(3) 

hp.write_map('/big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits', m=mm, coord='G', overwrite=True)

if not os.path.exists('/sptgrid/user/'+os.environ.get('USER')+'/'+pargs.group+'/'+mapname+'.fits'):
    os.system('gfal-copy -f /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/'+os.environ.get('USER')+'/'+pargs.group+'')
    os.system('rm /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits')

#mock-observe each map, coadd and calculate the power spectra
if pargs.submit:
#    if pargs.cleanup:
    if pargs.bkmask:
        os.system('python '+os.environ['HOME']+'/src/submit_'+pargs.group+'_2.1.py --inmap /sptgrid/user/'+os.environ['USER']+'/'+pargs.group+'/'+mapname+'.fits --submit --cleanup --bkmask --band '+pargs.band+' --nobs '+pargs.nobs+' --binw '+pargs.binw)   
    else:
        os.system('python '+os.environ['HOME']+'/src/submit_'+pargs.group+'_2.1.py --inmap /sptgrid/user/'+os.environ['USER']+'/'+pargs.group+'/'+mapname+'.fits --submit --cleanup --band '+pargs.band+' --nobs '+pargs.nobs+' --binw '+pargs.binw)
#    else:
#        os.system('python '+os.environ['HOME']+'/src/submit_'+pargs.group+'_2.1.py --inmap /sptgrid/user/'+os.environ['USER']+'/'+pargs.group+'/'+mapname+'.fits --submit --band '+pargs.band+' --nobs '+pargs.nobs)
else:
    os.system('python '+os.environ.get('HOME')+'/src/submit_'+pargs.group+'_2.1.py --inmap /sptgrid/user/'+os.environ['USER']+'/'+pargs.group+'/'+mapname+'.fits --band '+pargs.band+' --nobs '+pargs.nobs+' --binw '+pargs.binw)
