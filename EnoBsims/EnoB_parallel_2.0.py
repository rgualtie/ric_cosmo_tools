#EnoB
import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import sys
parser = argparse.ArgumentParser()
parser.add_argument('--nside', type=int, default=512,
                   help = 'Nside of the output maps, default=512')
parser.add_argument('--infile', type=str, default='/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat',
                   help = 'Filename of the .npy file containing the input spectra')
parser.add_argument('--band', type=str, default='90',
                   help = 'Frequency band the mock-observation, default=90')
parser.add_argument('--band2', type=str, default='150',
                   help = 'Frequency band the second mock-observation, default=1500')
parser.add_argument('--simidx', type=int, default=0,
                   help = 'Sim index, default=0')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of files used for the mock-observation, default=256')
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Delete the intermediate products')
pargs = parser.parse_args()

# example call: $python src/EnoB_parallel.py --submit --band 90 --band2 150 --simidx 99

if os.path.exists('/sptgrid/user/rgualtie/EnoB/output_spectra/spectra_EnoBlens_s1'+str(pargs.simidx).zfill(3)+'_idx'+str(pargs.simidx).zfill(3)+'_'+pargs.band+'GHz_cross_EnoBlens_s1'+str(pargs.simidx).zfill(3)+'_idx'+str(pargs.simidx).zfill(3)+'_'+pargs.band2+'GHz.npz'):
    print('output_spectrum for %i index already created, moving on' %pargs.simidx)
    sys.exit() 

f = pargs.infile
# make the maps
dls = np.genfromtxt(pargs.infile, unpack=True)[1:]
ell = np.arange(2,len(dls[0])+2)
lfac = ell*(ell+1)/2./np.pi
cls = dls/lfac

sd = 1000+pargs.simidx
np.random.seed(sd)
mm = hp.synfast(cls, nside=4*pargs.nside, lmax=2500, pol=True, new=True)

mapname = 'EnoBlens_s'+str(sd).zfill(4)+'_idx'+str(pargs.simidx).zfill(3) 

hp.write_map('/big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits', m=mm, coord='G', overwrite=True)

os.system('gfal-copy -f /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/'+os.environ.get('USER')+'/EnoB')
os.system('rm /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits')

#mock-observe each map, coadd and calculate the power spectra
if pargs.submit:
    if pargs.cleanup:
        os.system('python /home/'+os.environ.get('USER')+'/src/submit_EnoB_2.0.py --inmap /sptgrid/user/rgualtie/EnoB/'+mapname+'.fits --submit --cleanup --band '+pargs.band+' --band2 '+pargs.band2+' --nobs '+pargs.nobs)
    else:
        os.system('python /home/'+os.environ.get('USER')+'/src/submit_EnoB_2.0.py --inmap /sptgrid/user/rgualtie/EnoB/'+mapname+'.fits --submit --band '+pargs.band+' --band2 '+pargs.band2+' --nobs '+pargs.nobs)
else:
    os.system('python /home/'+os.environ.get('USER')+'/src/submit_EnoB_2.0.py --inmap /sptgrid/user/rgualtie/EnoB/'+mapname+'.fits --band '+pargs.band+' --band2 '+pargs.band2+' --nobs '+pargs.nobs)
