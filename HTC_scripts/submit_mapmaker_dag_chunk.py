import sys,os
import numpy as np
from glob import glob
import argparse, random
from spt3g import core, std_processing
from spt3g.cluster.condor_tools import condor_submit
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--dag', action = 'store_true',
                   help = 'Pack the jobs in a dag file')
parser.add_argument('--healpix', action = 'store_true',
                   help = 'Create and submit mock obseration healpix maps')
parser.add_argument('--nside',type=str, default=2048,
                   help = 'Healpix Nside for the output map')
parser.add_argument('-b','--band', action ='store',
                    default=None)
parser.add_argument('--chunkid',type=int, default=0,
                   help = 'Selects the data split, default = 0')
parser.add_argument('--chunklen',type=int, default=100,
                   help = 'Defines the chunk lenght, default = 100')
pargs = parser.parse_args()

# Example call:
# python submit_mapmaker_dag_chunk.py --submit --healpix --nside 2048 --band 90 --chunkid 0

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################
########## JOB AND USER SPECIFIC SETTINGS ########

# Description of job that will be appended to outputs
job_root = 'SPT3G_2019'
# Location of script to be run
script = os.path.join(spt3g_software,'std_processing/mapmakers/master_field_mapmaker.py')

grid_path = "/sptgrid/user/rgualtie/"
data_path = "/sptgrid/data/bolodata/downsampled/"

cal_files = sorted(glob("/sptgrid/analysis/calibration/calframe/ra0hdec-*/*.g3"))[pargs.chunkid*pargs.chunklen: (pargs.chunkid+1)*pargs.chunklen]
cal_files_no_path = []
for i in cal_files:
    cal_files_no_path.append(i.split('/')[-1])

input_files = sorted(glob(data_path+"/ra0hdec-*/*/[!nominal_online_cal]*"))[pargs.chunkid*pargs.chunklen: (pargs.chunkid+1)*pargs.chunklen]
obsids = []
field = []
for i in input_files:
    obsids.append(i.split('/')[-2])
    field.append(i.split('/')[-3])

aux_input_files = ['config_template_p10_%s.yaml'%pargs.band]

# Location for logging and data outputs
condor_dir = os.path.join('/scratch/rgualtie/condor_logs', job_root)
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

f = open(dag_script,'w')
for i,s in enumerate(input_files):
    jobname = job_root+'_'+field[i]+'_'+obsids[i]+'_'+s.split('/')[-1]
    if pargs.healpix:
        output_file = [jobname+'_hp.g3']
        extra_args = " --bands-to-use "+pargs.band+" -v --produce-simstub --produce-simstub --healpix --nside "+pargs.nside
    else:
        output_file = [jobname+'.g3']
        extra_args = " --bands-to-use "+pargs.band+" -v --produce-simstub --produce-simstub"
    try:
        idx = np.argwhere(np.array(cal_files_no_path) == s.split('/')[-2]+'.g3')[0][0]
    except:
        continue
    in_files = [s,cal_files[idx]]
    args = cal_files_no_path[idx]+" "+s.split('/')[-1]+" -o "+" ".join(output_file)+" --config-file "+" ".join(aux_input_files)+extra_args
    create_only = False
    if pargs.dag:
        create_only = True
    condor_submit(script, create_only=create_only, args = [args],
                  log_root = os.path.join(
                      condor_dir, field[i], obsids[i], str(pargs.band)+'GHz'),
                  output_root = os.path.join(
                      out_root, field[i], obsids[i], str(pargs.band)+'GHz'),
                  verbose = False,
                  retry = True,
                  jobname = jobname,
                  input_files = in_files,
                  aux_input_files=aux_input_files,
                  output_files=output_file,
                  requirements = requirements,
                  request_disk=request_disk,
                  request_memory=request_memory,
                  grid_proxy=grid_proxy)
    f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, field[i], obsids[i], str(pargs.band)+'GHz', jobname+'.submit')))
    f.write('RETRY %s 5\n'%jobname)
f.close()

if pargs.submit:
    os.system('condor_submit_dag -update_submit %s'%dag_script)
