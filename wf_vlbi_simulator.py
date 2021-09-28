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
commands.append('%s simulations/make_measurement_set.py single simulator_inputs.txt'%(inputs['stimela_exec']))

## Add noise to measurement sets & flag
commands.append('%s simulations/add_noise_hetero.py %s/single_pointing.ms %d'%(inputs['CASA_exec'],inputs['output_path'],int(inputs['size'])))

## Generate a terms
commands.append('%s simulations/generate_pb_aterms.py %s/single_pointing.ms 0 0 0'%(inputs['CASA_exec'],inputs['output_path']))

## Unzip a terms
commands.append('gunzip %s/single_pointing.ms_pb_flat_norotate.fits.gz'%(inputs['output_path']))

## Wsclean primary beam
commands.append('%s -name %s/single_pointing -no-update-model-required --aterm-kernel-size 157 -weight %s -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config single_pointing.ms_aterm_norotate_config.txt -size %d %d %s/single_pointing.ms'%(inputs['wsclean_exec'],inputs['output_path'],inputs['weight'],int(inputs['size']),int(inputs['size']),inputs['output_path']))

if inputs['mosaic'] == "False":
	commands.append('%s single_pointing-image-pb.fits'%(inputs['rms_exec']))
else:
	commands.append('%s simulations/generate_mosaic_pointings.py'%(inputs['CASA_exec']))

with open('job_%s.%s'%(step,params['job_manager']), 'a') as filehandle:
	for listitem in commands:
		filehandle.write('%s\n' % listitem)

if inputs['mosaic'] == "True":
	commands = []
	step = 'mosaic'
	write_hpc_headers(step,params)
