import os, hashlib
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--inmap',type=str, default='/sptgrid/user/rgualtie/BnoE/BnoElens_s1000_idx000.fits',
                   help = 'input map fits file, full path required')
parser.add_argument('--nside', type=str, default='512')
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--band2', type=str, default='150')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=256')
parser.add_argument('--cleanup', action = 'store_true',
                   help = 'Remove all the intermediate products')
pargs = parser.parse_args()

dag_script = '/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'.dag'

# Create the dag file
# Mock-obs the input map
os.system('python src/submit_mock_obs_map_dag.py --group BnoE --dag --inmap '+pargs.inmap+' --healpix --nside '+pargs.nside+' --band '+pargs.band+' --randomobs --nobs '+pargs.nobs+' --remove_duplicates')
os.system('python src/submit_mock_obs_map_dag.py --group BnoE --dag --inmap '+pargs.inmap+' --healpix --nside '+pargs.nside+' --band '+pargs.band2+' --randomobs --nobs '+pargs.nobs+' --remove_duplicates')

# Coadd the fields
coadd_mapname = os.path.splitext(os.path.basename(pargs.inmap))[0]

os.system('python src/submit_coadder_stdproc.py --nobs '+pargs.nobs+' --inmap '+coadd_mapname+' --dag --group BnoE --band '+pargs.band)
os.system('python src/submit_coadder_stdproc.py --nobs '+pargs.nobs+' --inmap '+coadd_mapname+' --dag --group BnoE --band '+pargs.band2)

os.system('python src/submit_cross_spectra.py --inmap /sptgrid/user/'+os.environ.get('USER')+'/BnoE/'+coadd_mapname+'/'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3 --inmap2 /sptgrid/user/'+os.environ.get('USER')+'/BnoE/'+coadd_mapname+'/'+coadd_mapname+'_'+pargs.band2+'GHz_nstubs'+pargs.nobs+'.g3 --dag --group BnoE --band '+pargs.band+' --band2 '+pargs.band2+' --nobs '+pargs.nobs)

# Define the dependencies
# mock-obs job names
jn = np.genfromtxt('/big_scratch/'+os.environ.get('USER')+'/condor_logs/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'/'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_job_names.txt', unpack=True, dtype=str)

with open(dag_script,'a') as f:
    f.write('PARENT '+' '.join(jn)+' CHILD coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'\n')
    f.write('PARENT '+' '.join(jn)+' CHILD coadd_'+coadd_mapname+'_'+pargs.band2+'GHz_nstubs'+pargs.nobs+'\n')
    f.write('PARENT coadd_'+coadd_mapname+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' coadd_'+coadd_mapname+'_'+pargs.band2+'GHz_nstubs'+pargs.nobs+' CHILD spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_cross_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band2+'GHz\n')
    if pargs.cleanup:
        f.write('SCRIPT POST spectra_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band+'GHz_cross_'+os.path.splitext(os.path.basename(pargs.inmap))[0]+'_'+pargs.band2+'GHz /home/'+os.environ.get('USER')+'/src/cleanup_BnoE.sh '+os.path.splitext(os.path.basename(pargs.inmap))[0])

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
