#mix matrix
import numpy as np
import os, sys, argparse
import healpy as hp
from glob import glob
from spt3g import core, mapmaker, util, mapspectra, maps, std_processing

parser = argparse.ArgumentParser()
parser.add_argument('--nside', type=int, default=512,
                   help = 'Nside of the output maps, default=512')
parser.add_argument('--infile', type=str, 
                   default='/sptgrid/user/rgualtie/mix_matrix/input_pulse_spectra/pulse_spectra_EE_0033.npy',
                   help = 'Filename of the .npy file containing the input spectra')
parser.add_argument('--band', type=str, default='90',
                   help = 'Frequency band the mock-observation, default=90')
parser.add_argument('--simidx', type=int, default=0,
                   help = 'Sim index, default=0')
parser.add_argument('--submit', action = 'store_true', 
                   help = 'Submit the dag file')
pargs = parser.parse_args()

# example call: $python src/mix_matrix_parallel_2.0.py --submit --nside 512 --infile /sptgrid/user/rgualtie/mix_matrix/input_spectra/spectra_TT_0008.npy --band 150 --simidx 99
for s in ['TT','EE','BB']:
    if s in pargs.infile:
        sp = s 
sbin = pargs.infile.split('/')[-1].replace('.npy','')[-4:]
if os.path.exists('/sptgrid/user/rgualtie/mix_matrix/output_pulse_spectra/spectra_'+sp+'_'+sbin+'_s1'+str(pargs.simidx).zfill(3)+'_idx'+str(pargs.simidx).zfill(3)+'_'+pargs.band+'GHz_nstubs256.npz'):
    print('output_spectrum for %s already created, moving on' %pargs.infile.split('/')[-1])
    sys.exit()
f = pargs.infile
# make the maps
if os.path.exists(pargs.infile):
    dls = np.load(pargs.infile)
else:
    os.system('python $HOME/src/make_mixmat_pulse_input_spectra.py')
    dls = np.load(pargs.infile)
sd = 1000+pargs.simidx
np.random.seed(sd)
ell = np.arange(len(dls[0]))
lfac = ell*(ell+1)/2./np.pi+1e-6
cls = dls/lfac
mm = hp.synfast(cls, nside=pargs.nside, pol=True, new=True)

mapname = os.path.splitext(os.path.basename(pargs.infile))[0]+'_s'+str(sd).zfill(4)+'_idx'+str(pargs.simidx).zfill(3) 
mapname = mapname.replace('pulse_spectra_','')
hp.write_map('/big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits', m=mm, coord='G', overwrite=True)
os.system('gfal-copy -f /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/'+os.environ.get('USER')+'/mix_matrix')
os.system('rm /big_scratch/'+os.environ.get('USER')+'/'+mapname+'.fits')

#mock-observe each map, coadd and calculate the power spectra
if pargs.submit:
    os.system('python /home/'+os.environ.get('USER')+'/src/submit_mixmatrix.py --inmap /sptgrid/user/rgualtie/mix_matrix/'+mapname+'.fits --submit --cleanup --band '+pargs.band)
else:
    os.system('python /home/'+os.environ.get('USER')+'/src/submit_mixmatrix.py --inmap /sptgrid/user/rgualtie/mix_matrix/'+mapname+'.fits --cleanup --band '+pargs.band)
