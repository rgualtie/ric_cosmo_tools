import numpy as np
import os, sys, time
import argparse
from glob import glob
from spt3g import core, std_processing
from spt3g.cluster.condor_tools import condor_submit

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--nobs', type=int, default=64,
                   help = 'Number of stubs used to generate the file, default = 64')
parser.add_argument('--band', type=str, default='90',
                   help = 'Band of the map, default = 90')
parser.add_argument('--dag', action = 'store_true',
                   help = 'Create jobs and also submit them through a dag file')
parser.add_argument('--inmap',type=str,
                   help = 'input maps path on the grid, full path is not required')
parser.add_argument('--group',type=str,default='',
                   help = 'Group the output to a common folder')
pargs = parser.parse_args()

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

mapname = pargs.inmap
# coadd the mock observed map
# where the mock observed g3 files live

grid_path = os.path.join('/sptgrid/user/rgualtie',pargs.group)
#print(grid_path)
script = os.path.join(spt3g_software,'std_processing/combining_maps.py') 
job_root = mapname

# Location for logging and data outputs
condor_dir = os.path.join('/big_scratch/rgualtie/condor_logs', job_root)
out_root = os.path.join('/sptgrid/user/rgualtie',pargs.group, job_root)

# Grid proxy location
grid_proxy = "/home/rgualtie/rgualtie_proxy"

# Hardware requirements
request_disk = 16*core.G3Units.GB
request_memory = 16*core.G3Units.GB

# These requirements were adapted from a similar string in doc/osg/osg_guide.md
# It prevents the same machine from being used on re-submissions,
# requires remote machines to use RHEL 7, and forbids NPX.
requirements = '''( ( HAS_CVMFS_spt_opensciencegrid_org ) && ( ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4 ) || ( RCC_Factory == "ciconnect" ) ) && ( ( TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName5 ) || ( RCC_Factory == "ciconnect" ) ) ) && ( OSGVO_OS_STRING == "RHEL 7" ) && (GLIDEIN_ResourceName =!= "NPX"))'''

# Create and submit jobs
if not os.path.exists(condor_dir):
    os.mkdir(condor_dir)

# What a file looks like:
# /sptgrid/user/rgualtie/EnoB/EnoBlens_s1031_idx031_90GHz_ra0hdec-52.25_90555480_hp_nstubs256.g3

input_files = glob(grid_path+'/'+mapname+'/ra0hdec-*/*/'+pargs.band+'GHz/*nstubs'+str(pargs.nobs)+'.g3')

if not input_files:
    input_files = []
    jns = np.genfromtxt(condor_dir+'/'+mapname+'_job_names.txt', unpack=True, dtype=str)
    # Files look like: EnoBlens_s1000_idx000_90GHz_ra0hdec-44.75_80887441
    for field in ['44.75','52.25','59.75','67.25']:
        obsids = []
        for jn in jns:
            if pargs.band+'GHz' in jn:
                if field in jn:
                    obsids.append(jn.split('_')[-1])
        for obsid in obsids:
                input_files.append(grid_path+'/'+mapname+'/ra0hdec-'+field+'/'+obsid+'/'+pargs.band+'GHz/'+mapname+'_'+pargs.band+'GHz_ra0hdec-'+field+'_'+obsid+'_hp_nstubs'+str(pargs.nobs)+'.g3')

base_input_files = [os.path.basename(infile) for infile in input_files]
jobname = 'coadd_'+mapname+'_'+pargs.band+'GHz_nstubs'+str(pargs.nobs)
args = " ".join(base_input_files)+" -o "+jobname+".g3"
create_only = False
if pargs.dag:
    create_only = True

condor_submit(script, create_only=create_only, args = [args],
          log_root = os.path.join(condor_dir),
          output_root = os.path.join(out_root),
          verbose = False,
          retry = False,
          jobname = jobname,
          input_files = input_files,
          ignore_input = True,
          aux_input_files=[],
          output_files = [jobname+'.g3'],
          requirements = requirements,
          request_disk=request_disk,
          request_memory=request_memory,
          grid_proxy=grid_proxy,
          force_build_tarball=False)

if pargs.dag:
    dag_script = os.path.join(condor_dir, job_root+'.dag')
    with open(dag_script,'a') as f:
        f.write('JOB %s %s\n'%(jobname, condor_dir+'/'+jobname+'.submit'))
        f.write('RETRY %s 5\n'%jobname)
    
    if pargs.submit:
        os.system('condor_submit_dag -update_submit %s'%dag_script)
