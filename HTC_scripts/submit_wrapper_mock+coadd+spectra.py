import os, hashlib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--inmap',type=str, default='/sptgrid/user/rgualtie/PlanckDR3Maps/HFI_SkyMap_100_2048_R3.01_full_C.fits',
                   help = 'input map fits file, full path required')
parser.add_argument('--nside', type=str, default='512')
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--nobs', type=str, default='128',
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=128')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Remove all the intermediate products')
pargs = parser.parse_args()

dag_script = '/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'.dag'

# Create the dag file
# Mock-obs the input map
os.system('python src/submit_mock_obs_map_dag.py --dag --inmap '+pargs.inmap+' --healpix --nside '+pargs.nside+' --band '+pargs.band+' --randomobs --nobs '+pargs.nobs)
# Coadd the fields
coadd_mapname = os.path.splitext(os.path.basename(pargs.inmap))[0]
os.system('python src/submit_coadder_stdproc.py --nobs '+pargs.nobs+' --inmap '+coadd_mapname+' --dag --band '+pargs.band)
os.system('python src/submit_power_spectra.py --inmap /sptgrid/user/'+os.environ.get('USER')+'/'+coadd_mapname+'/'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3 --dag')

# Define the dependencies
# mock-obs job names
jn = np.genfromtxt('/scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_job_names.txt', unpack=True, dtype=str)

with open(dag_script,'a') as f:
    f.write('PARENT '+' '.join(jn)+' CHILD coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs)
    f.write('\nPARENT coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' CHILD spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz\n')
    if pargs.cleanup:
        #os.system('gfal-copy /home/'+os.environ.get('USER')+'/src/cleanup.sh gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/'+os.environ.get('USER'))
        f.write('SCRIPT POST spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz /home/'+os.environ.get('USER')+'/src/cleanup.sh '+os.path.splitext(os.path.basename(pargs.inmap))[0])

# Cleanup dag file duplicate lines
output_file_path = '/scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_clean.dag'
output_file = open(output_file_path, "w")
completed_lines_hash = set()
for line in open(dag_script, "r"):
    hashValue = hashlib.md5(line.rstrip().encode('utf-8')).hexdigest()
    if hashValue not in completed_lines_hash:
        output_file.write(line)
        completed_lines_hash.add(hashValue)
output_file.close()
    
# Submit the dag
if pargs.submit:
    os.system('condor_submit_dag -update_submit %s'%output_file_path)
