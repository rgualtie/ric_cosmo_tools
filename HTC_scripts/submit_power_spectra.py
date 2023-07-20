import numpy as np
import os, sys, time
import argparse
from glob import glob
from spt3g import core, std_processing
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--dag', action = 'store_true',
                   help = 'Create jobs and also submit them through a dag file')
parser.add_argument('--inmap',type=str,
                    help = 'input maps path on the grid, full path is required')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of stubs used to generate the file, default = 256')
parser.add_argument('--band', type=str, default='90',
                   help = 'Band of the map, default = 90')
parser.add_argument('--group',type=str, default='',
                   help = 'Group the output in the selected folder')
pargs = parser.parse_args()
try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

# where the coadd g3 file live
grid_path = os.path.join('/sptgrid/user/rgualtie', pargs.group)

# spectra estimator wrapper
script = '/home/rgualtie/src/make_spectra.py' 
job_root = os.path.dirname(pargs.inmap.replace('coadd_',''))

# Location for logging and data outputs
# clean the path name to match the input map basename with no extensions
map_name = os.path.basename(pargs.inmap).replace('coadd_','').replace('_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.g3','')

condor_dir = os.path.join('/big_scratch/rgualtie/condor_logs', map_name)

out_root = os.path.join(grid_path, job_root)

dag_script = os.path.join(condor_dir, map_name+'.dag')

# Grid proxy location
grid_proxy = "/home/rgualtie/rgualtie_proxy"

# Hardware requirements
request_disk = 8*core.G3Units.GB
request_memory = 8*core.G3Units.GB

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

# Create and submit jobs
if not os.path.exists(condor_dir):
    os.mkdir(condor_dir)

if pargs.dag:
    f = open(dag_script,'a')

base_input_file = os.path.basename(pargs.inmap) 

mask = '/sptgrid/user/rgualtie/ApodMaskPS2048.fits'
jobname = 'spectra_'+map_name+'_'+pargs.band+'GHz_nstubs'+pargs.nobs
args = "--inmap coadd_"+base_input_file+" --mask "+os.path.basename(mask)+" --output "+jobname+'.npz'

create_only = True 
if pargs.submit:
    create_only = False
if pargs.dag:
    create_only = True

condor_submit(script, create_only=create_only, args = [args],
          log_root = os.path.join(condor_dir),
          output_root = os.path.join(grid_path, 'output_spectra/'+pargs.band+'GHz'),
          verbose = False,
          retry = True,
          jobname = jobname,
          input_files = [job_root+'/coadd_'+base_input_file],
          ignore_input = True,
          aux_input_files=[mask],
          output_files = [jobname+'.npz'],
          requirements = requirements,
          request_disk=request_disk,
          request_memory=request_memory,
          grid_proxy=grid_proxy)

if pargs.dag:
    f.write('JOB %s %s\n'%(jobname, condor_dir+'/'+jobname+'.submit'))
    f.write('RETRY %s 5\n'%jobname)
    f.close()
    if pargs.submit:
        os.system('condor_submit_dag -update_submit %s'%dag_script)
