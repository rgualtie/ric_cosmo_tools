#mix matrix
import numpy as np
import os, time, argparse
import healpy as hp
from glob import glob
import subprocess
from spt3g import core, mapmaker, util, mapspectra, maps, std_processing

parser = argparse.ArgumentParser()
parser.add_argument('--nside', type=int, default=512,
                   help = 'Nside of the output maps, default=512')
parser.add_argument('--infile', type=str, 
                   default='/sptgrid/user/rgualtie/mix_matrix/input_spectra/spectra_EE_0033.npy',
                   help = 'Filename of the .npy file containing the input spectra')
parser.add_argument('--nobs', type=int, default=64,
                   help = 'Number of ObsIds for the mock-observation, default=64')
parser.add_argument('--simidx', type=int, default=0,
                   help = 'Sim index, default=0')
parser.add_argument('--debug', action = 'store_true',
                   help = 'Turns on the debug mode, it inhibits the cleanup of the intermediate products')
pargs = parser.parse_args()

f = pargs.infile
# make the maps
if os.path.exists(pargs.infile):
    cls = np.load(pargs.infile)
else:
    os.system('python $HOME/src/make_mixmat_input_spectra.py')
    cls = np.load(pargs.infile)
sd = 1000+pargs.simidx
np.random.seed(sd)
mm = hp.synfast(cls, nside=pargs.nside, pol=True, new=True)
mapname = f.split('/')[-1][:-4]+'_s'+str(sd).zfill(4)+'_idx'+str(pargs.simidx).zfill(3)
hp.write_map(mapname+'.fits', m=mm, coord='G', overwrite=True)
os.system('gfal-copy '+mapname+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix')
os.system('rm '+mapname+'.fits')
#mock-observe each map
os.system('for f in 90 150 220;do python $HOME/src/submit_mock_obs_map_dag.py --submit --inmap /sptgrid/user/rgualtie/mix_matrix/'+mapname+'.fits --healpix --nside '+str(pargs.nside)+' --band $f --randomobs --nobs '+str(pargs.nobs)+' --mixmat;done')
# Wait until the mock is done
while 1:
    os.system('condor_q -nobatch > temp.txt')
    queue = np.genfromtxt('temp.txt', unpack=True, dtype=str, delimiter='\t') 
    queue = queue[2:-3]
    running_procs = []
    for i in range(len(queue)):
        running_procs.append(queue[i].split()[-1])
    if any(mapname in s for s in running_procs):
        # wait 10min and check again
        print('The grid is still processing the map..check in 10min')
        os.system('rm temp.txt')
        time.sleep(600)
    else:
        print('%s is ready for coadd'%mapname)
        break
# Coadd the mock observed maps, generate the spectra and cleanup
if pargs.debug:
    os.system('python $HOME/src/coadder.py --inmap '+mapname+' --nobs '+str(pargs.nobs)+' --chunk_len 64 --spectra --mixmat')
else:
    os.system('python $HOME/src/coadder.py --inmap '+mapname+' --nobs '+str(pargs.nobs)+' --chunk_len 64 --spectra --mixmat --cleanup')
    os.system('gfal-rm gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/'+mapname+'.fits')
