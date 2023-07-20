import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--inmap',type=str, default='HFI_SkyMap_100_2048_R3.01_full_C',
                   help = 'input map fits file, full path required')
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--nobs', type=str, default='128',
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=128')
pargs = parser.parse_args()

dag_script = '/big_scratch/rgualtie/condor_logs/'+pargs.inmap+'/'+pargs.inmap+'.dag'

# Create the dag file
os.system('python src/submit_coadder_stdproc.py --nobs '+pargs.nobs+' --inmap '+pargs.inmap+' --dag')
os.system('python src/submit_power_spectra.py --inmap /sptgrid/user/rgualtie/'+pargs.inmap+'/'+pargs.inmap+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3 --dag')

# Define the dependencies
with open(dag_script,'a') as f:
    f.write('PARENT coadd_'+pargs.inmap+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+' CHILD spectra_coadd_'+pargs.inmap+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz')

# Submit the dag
if pargs.submit:
    os.system('condor_submit_dag -update_submit %s'%dag_script)
