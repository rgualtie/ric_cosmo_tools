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
                    help = 'input maps path, full path is not required')
parser.add_argument('--nobs', type=int, default=64,
                   help = 'Number of stubs to be used')
pargs = parser.parse_args()

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

mapname = pargs.inmap
# coadd the mock observed map
# where the mock observed g3 files live
grid_path = '/sptgrid/user/rgualtie/'

script = '/home/rgualtie/src/coadder.py'
job_root = 'coadd_'+mapname

# Location for logging and data outputs
condor_dir = os.path.join('/big_scratch/rgualtie/condor_logs', job_root)

out_root = os.path.join('/sptgrid/user/rgualtie', job_root)

dag_script = os.path.join(condor_dir, job_root+ '.dag')

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
    f = open(dag_script,'w')

#input_files =  
args = '--inmap '+pargs.inmap+' --nobs '+str(pargs.nobs)+' --chunk_len 64 --spectra --cleanup' 
jobname = 'coadd_'+mapname

create_only = False
if pargs.dag:
    create_only = True
condor_submit(script, create_only=create_only, args = [args],
              log_root = os.path.join(condor_dir),
              output_root = os.path.join(out_root),
              verbose = False,
              retry = False,
              jobname = jobname,
              #input_files=input_files,
              output_files=[jobname+'.g3'],
              requirements = requirements,
              request_disk=request_disk,
              request_memory=request_memory,
              grid_proxy=grid_proxy,
              force_build_tarball=False)
if pargs.dag:
    f.write('JOB %s %s\n'%(jobname, jobname+'.submit'))
    f.write('RETRY %s 1\n'%jobname)

if pargs.dag:
    f.close()
    if pargs.submit:
        os.system('condor_submit_dag -update_submit -maxidle 25 %s'%dag_script)
