import sys,os
import numpy as np
from glob import glob
import argparse, random, time
from spt3g import core#, std_processing
from spt3g.cluster.condor_tools import condor_submit
from datetime import datetime
import hashlib 

parser = argparse.ArgumentParser()
parser.add_argument('--submit', action = 'store_true',
                   help = 'Create jobs and also submit them')
parser.add_argument('--dag', action = 'store_true',
                   help = 'Create jobs and also submit them through a dag file')
parser.add_argument('--remove_duplicates', action = 'store_true',
                   help = 'Cleanup the dag file from duplicate lines')
parser.add_argument('--inmap',type=str, 
                   help = 'input map fits file, full path required')
parser.add_argument('--healpix', action = 'store_true',
                   help = 'Create and submit mock obseration healpix maps')
parser.add_argument('--nside',type=str, default=2048,
                   help = 'Healpix Nside for the output map')
parser.add_argument('-b','--band', action ='store',
                    default=None)
parser.add_argument('--randomobs', action ='store_true',
                   help = 'If present selects a subset of ObsIDs, number is defined by --nobs 64')
parser.add_argument('--nobs', type=int, default=256,
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=256')
parser.add_argument('--group',type=str, default='',
                   help = 'Group the output into a single folder')
pargs = parser.parse_args()

try:
    spt3g_software = os.environ['SPT3G_SOFTWARE_PATH']
except KeyError as e:
    print('%s environment variable not set'%e)
    raise

##################################################
########## JOB AND USER SPECIFIC SETTINGS ########

# Description of job that will be appended to outputs
map_name = pargs.inmap.split('/')[-1]
job_root = map_name.replace('.fits','') 
#current date and time
now = datetime.now()
format = "%d%m%Y_%H%M_"
#format datetime using strftime()
time1 = now.strftime(format)
# Location of script to be run
script = os.path.join(spt3g_software,'std_processing/mapmakers/master_field_mapmaker.py')

stubs = glob('/sptgrid/user/javva/full_2019_and_2020/ra0hdec-*/full_2019_p10_elnod_pair_difference_healpix_recent_cal_lpf3k/*/simstub_*.g3')
stubs44 = glob('/sptgrid/user/javva/full_2019_and_2020/ra0hdec-44.75/full_2019_p10_elnod_pair_difference_healpix_recent_cal_lpf3k/*/simstub_*.g3')
stubs52 = glob('/sptgrid/user/javva/full_2019_and_2020/ra0hdec-52.25/full_2019_p10_elnod_pair_difference_healpix_recent_cal_lpf3k/*/simstub_*.g3')
stubs59 = glob('/sptgrid/user/javva/full_2019_and_2020/ra0hdec-59.75/full_2019_p10_elnod_pair_difference_healpix_recent_cal_lpf3k/*/simstub_*.g3')
stubs67 = glob('/sptgrid/user/javva/full_2019_and_2020/ra0hdec-67.25/full_2019_p10_elnod_pair_difference_healpix_recent_cal_lpf3k/*/simstub_*.g3')

if pargs.randomobs:
    random.seed(random.randint(1,100)*int(time.time()))
    stubs44 = random.choices(stubs44, k=pargs.nobs//4)
    stubs52 = random.choices(stubs52, k=pargs.nobs//4)
    stubs59 = random.choices(stubs59, k=pargs.nobs//4)
    stubs67 = random.choices(stubs67, k=pargs.nobs//4)
    #stubs = random.choices(stubs, k=pargs.nobs)
    stubs = stubs44+stubs52+stubs59+stubs67
 
no_path_stubs = []
field = []
obsids = []
for i in stubs: 
    no_path_stubs.append(i.split('/')[-1])
    field.append(i.split('/')[-4])
    obsids.append(i.split('/')[-2])

aux_input_files = ['config_template_p10_healpix_lpf3k_%sGHz.yaml'%pargs.band]
# Command line args to be passed to the script on remote machine.

# Location for logging and data outputs
condor_dir = os.path.join('/big_scratch/'+os.environ.get('USER')+'/condor_logs', job_root)

out_root = os.path.join('/sptgrid/user/'+os.environ.get('USER'), pargs.group, job_root)
#print(out_root)
dag_script = os.path.join(condor_dir, job_root+ '.dag')

# Grid proxy location
grid_proxy = '/home/'+os.environ.get('USER')+'/'+os.environ.get('USER')+'_proxy'

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

#if os.path.exists(condor_dir+'/'+job_root+'_job_names.txt'):
#    os.system('rm '+condor_dir+'/'+job_root+'_job_names.txt')

if pargs.dag:
    f = open(dag_script,'a')
for i,s in enumerate(stubs):
    jobname = job_root+'_'+pargs.band+'GHz_'+field[i]+'_'+obsids[i]
    with open(condor_dir+'/'+job_root+'_job_names.txt','a') as jn:
        jn.write(jobname+'\n')
    if pargs.healpix:
        if pargs.nobs:
            output_file = [jobname+'_hp_nstubs'+str(pargs.nobs)+'.g3']
        else:
            output_file = [jobname+'_hp.g3']
        extra_args = "--bands-to-use "+pargs.band+" --healpix --nside "+pargs.nside
    
    else:
        output_file = [jobname+'.g3']
        extra_args = "--bands-to-use "+pargs.band
    input_files = [s, pargs.inmap]
    args = no_path_stubs[i]+" -o "+" ".join(output_file)+" --config-file "+" ".join(aux_input_files)+" --sim --sim-map "+map_name+" "+extra_args
    sky_sim_file = map_name
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
                  input_files=input_files,
                  ignore_input=True,
                  aux_input_files=aux_input_files,
                  output_files=output_file,
                  requirements = requirements,
                  request_disk=request_disk,
                  request_memory=request_memory,
                  grid_proxy=grid_proxy,
                  force_build_tarball=False)
    if pargs.dag:
        f.write('JOB %s %s\n'%(jobname, os.path.join(condor_dir, field[i], obsids[i], str(pargs.band)+'GHz', jobname+'.submit')))
        f.write('RETRY %s 50\n'%jobname)

# Cleanup dag file duplicate lines
if pargs.remove_duplicates:
    output_file_path = dag_script[:-4]+'_clean.dag'
    clean_dag = open(output_file_path, "w")
    completed_lines_hash = set()
    for line in open(dag_script, "r"):
        hashValue = hashlib.md5(line.rstrip().encode('utf-8')).hexdigest()
        if hashValue not in completed_lines_hash:
            clean_dag.write(line)
            completed_lines_hash.add(hashValue)
    clean_dag.close()
if pargs.dag:
    f.close()
    if pargs.submit:
        if pargs.remove_duplicates:
            os.system('condor_submit_dag -maxidle 50 -update_submit %s'%output_file_path)
        else:
            os.system('condor_submit_dag -maxidle 50 -update_submit %s'%dag_script)
