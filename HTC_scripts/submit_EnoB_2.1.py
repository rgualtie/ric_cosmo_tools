import os, hashlib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--inmap',type=str, default='/sptgrid/user/'+os.environ['USER']+'/EnoB/EnoBlens_s1000_idx000.fits',
                   help = 'input map fits file, full path required')
parser.add_argument('--nside', type=str, default='2048')
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=256')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Remove all the intermediate products')
parser.add_argument('--bkmask', action = 'store_true',
                   help = 'Uses bk14 mask')
parser.add_argument('--binw', type=str, default='25')
pargs = parser.parse_args()

dag_script = '/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz.dag'

# Create the dag file
# Mock-obs the input map
os.system('python '+os.environ['HOME']+'/src/submit_mock_obs_map_dag_2.0.py --group EnoB --dag --inmap '+pargs.inmap+' --healpix --nside '+pargs.nside+' --band '+pargs.band+' --randomobs --nobs '+pargs.nobs+' --remove_duplicates')

# Coadd the fields
coadd_mapname = os.path.splitext(os.path.basename(pargs.inmap))[0]

os.system('python '+os.environ['HOME']+'/src/submit_coadder_stdproc_2.0.py --nobs '+pargs.nobs+' --inmap '+coadd_mapname+' --dag --group EnoB --band '+pargs.band)

if pargs.bkmask:
    os.system('python '+os.environ['HOME']+'/src/submit_power_spectra_2.0.py --inmap /sptgrid/user/'+os.environ.get('USER')+'/EnoB/'+coadd_mapname+'/'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3 --dag --group EnoB --bkmask --band '+pargs.band+' --nobs '+pargs.nobs+' --binw '+pargs.binw)
else:
    os.system('python '+os.environ['HOME']+'/src/submit_power_spectra_2.0.py --inmap /sptgrid/user/'+os.environ.get('USER')+'/EnoB/'+coadd_mapname+'/'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3 --dag --group EnoB --band '+pargs.band+' --nobs '+pargs.nobs+' --binw '+pargs.binw)

# Define the dependencies
# mock-obs job names
jn = np.genfromtxt('/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_job_names.txt', unpack=True, dtype=str)

with open(dag_script,'a') as f:
    f.write('PARENT '+' '.join(jn)+' CHILD coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs)
    f.write('\nPARENT coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' CHILD spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'\n')
    if pargs.cleanup:
        f.write('SCRIPT POST spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' '+os.environ.get('HOME')+'/src/cleanup_EnoB_2.0.sh '+os.path.splitext(os.path.basename(pargs.inmap))[0]+' '+pargs.band)

# Cleanup dag file duplicate lines
#output_file_path = dag_script[:-4]+'_clean.dag'
output_file_path = '/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_clean.dag'
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
