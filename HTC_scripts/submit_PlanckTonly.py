import os, hashlib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--inmap',type=str, default='/sptgrid/user/rgualtie/PlanckTonlyMaps/P100Tonly_full.fits',
                   help = 'input map fits file, full path required')
parser.add_argument('--nside', type=str, default='512')
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--nobs', type=str, default='512',
                   help = 'Number of ObsIDs to be used, default=512')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Remove all the intermediate products')
pargs = parser.parse_args()

# example call: $python ~/src/submit_PlanckTonly.py --submit --inmap /sptgrid/user/rgualtie/PlanckTonlyMaps/P100Tonly_full.fits --nside 512 --band 90 --nobs 512 --cleanup

dag_script = '/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'.dag'

# Create the dag file
# Mock-obs the input map
os.system('python src/submit_mock_obs_map_dag.py --group PlanckTonlyMaps --dag --inmap '+pargs.inmap+' --healpix --nside '+pargs.nside+' --band '+pargs.band+' --randomobs --nobs '+pargs.nobs)

# Coadd the fields
coadd_mapname = os.path.splitext(os.path.basename(pargs.inmap))[0]

os.system('python src/submit_coadder_stdproc.py --nobs '+pargs.nobs+' --inmap '+coadd_mapname+' --dag --group PlanckTonlyMaps --band '+pargs.band)

# Define the dependencies
# mock-obs job names
jn = np.genfromtxt('/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_job_names.txt', unpack=True, dtype=str)

with open(dag_script,'a') as f:
    f.write('PARENT '+' '.join(jn)+' CHILD coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'\n')
    if pargs.cleanup:
        f.write('SCRIPT POST coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' /home/'+os.environ.get('USER')+'/src/cleanup_Tonly.sh '+os.path.splitext(os.path.basename(pargs.inmap))[0]+' '+pargs.band+' '+pargs.nobs)

# Cleanup dag file duplicate lines
output_file_path = dag_script[:-4]+'_clean.dag'
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
    #os.system('condor_submit_dag -update_submit %s'%dag_script)
