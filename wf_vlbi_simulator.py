import inspect, os, sys, json, ast
import copy
## Python 2 will need to adjust for casa 6
import collections

filename = inspect.getframeinfo(inspect.currentframe()).filename
sys.path.append(os.path.dirname(os.path.realpath(filename)))
sys.path.append(os.path.dirname(os.path.realpath(filename))+"/simulations")

from simulator_functions import *
from wf_vlbi_functions import headless

## Imports input_file
try:
	i = sys.argv.index("-c") + 2
except:
	i = 1
	pass

## Load global inputs
inputs = headless(sys.argv[i])

## Set the parameters for the HPC resources
params = {}
params['job_manager'] = str(inputs['job_manager'])
params['email_progress'] = str(inputs['email_progress'])
params['HPC_project_code'] = str(inputs['HPC_project_code'])
params['partition'] = str(inputs['partition'])
params['walltime'] = str(inputs['walltime'])
params['nodetype'] = str(inputs['nodetype'])
params['nodes'] = int(inputs['nodes'])
params['cpus'] = int(inputs['cpus'])
params['mpiprocs'] = int(inputs['mpiprocs'])

## Set the parameters for software
print(inputs['CASA_exec'])
'''
singularity exec /idia/software/containers/casa-6.3.simg python make_itrf.py evn emerlin
singularity exec ../STIMELA_SINGULARITY_IMAGES/stimela_simms_1.2.0.img python make_measurement_set.py single evn emerlin
singularity exec /idia/software/containers/casa-6.3.simg python add_noise_hetero.py evn_emerlin.ms 2048
singularity exec /idia/software/containers/casa-6.3.simg python generate_pb_aterms.py evn_emerlin.ms 0 0 0 
gunzip "evn_emerlin.ms_pb_flat_norotate.fits.gz"
singularity exec /idia/software/containers/wsclean-gpu.simg wsclean -name evn_emerlin_nat -no-update-model-required --aterm-kernel-size 157 -weight natural -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config evn_emerlin.ms_aterm_norotate_config.txt -size 2048 2048 evn_emerlin.ms
singularity exec /idia/software/containers/wsclean-gpu.simg wsclean -name evn_emerlin_uni -no-update-model-required --aterm-kernel-size 157 -weight uniform -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config evn_emerlin.ms_aterm_norotate_config.txt -size 2048 2048 evn_emerlin.ms
'''
## Generate single pointing to fit beam
commands = []
step = 'single_pointing'
write_hpc_headers(step,params)

## Generate itrfs
antennae = ast.literal_eval(inputs['antennae'])
commands.append('%s simulations/make_itrf.py %s'%(inputs['CASA_exec']," ".join(antennae)))

## Generate measurement set
commands.append('%s simulations/make_measurement_set.py single %s %s'%(inputs['stimela_exec'],inputs['output_path']," ".join(antennae)))

with open('job_%s.%s'%(step,params['job_manager']), 'a') as filehandle:
	for listitem in commands:
		filehandle.write('%s\n' % listitem)
