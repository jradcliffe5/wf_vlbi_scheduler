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
part = int(sys.argv[i+1])

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
params['output_path'] = inputs['output_path']
params['mem'] = inputs['mem']
params['max_jobs'] = int(inputs['max_jobs'])

## Set the parameters for software
print(inputs['CASA_exec'])

obs_freq = float(inputs['obs_freq'])

if (obs_freq > 1.0) & (obs_freq < 2.0):
	band='L'
elif (obs_freq > 2.0) & (obs_freq < 3.5):
	band='S'
elif (obs_freq > 3.5) & (obs_freq < 7.5):
	band='C'
elif (obs_freq > 7.5):
	band='K'
else:
	print('band not supported')
	sys.exit()

## Generate single pointing to fit beam
if part == 1:
	commands = []
	step = 'single_pointing'
	write_hpc_headers(step,params)
	#commands.append('module purge singularity')

	## Generate itrfs
	antennae = ast.literal_eval(inputs['antennae'])
	commands.append('%s simulations/make_itrf.py %s'%(inputs['CASA_exec']," ".join(antennae)))

	## Generate measurement set
	commands.append('%s simulations/make_measurement_set.py single simulator_inputs.txt'%(inputs['stimela_exec']))

	## Add noise to measurement sets & flag
	commands.append('%s simulations/add_noise_hetero.py %s/single_pointing.ms %d %.3f %s %s'%(inputs['CASA_exec'],inputs['output_path'],int(inputs['size']),float(inputs['time_multiplier']),band,inputs['cell']))

	## Generate a terms
	commands.append('%s simulations/generate_pb_aterms.py %s/single_pointing.ms 0 0 0 %s'%(inputs['CASA_exec'],inputs['output_path'],band))

	## Unzip a terms
	commands.append('gunzip -f %s/single_pointing.ms_pb_flat_norotate.fits.gz'%(inputs['output_path']))

	## Wsclean primary beam
	commands.append('%s -name %s/single_pointing -no-update-model-required --aterm-kernel-size 157 -weight %s -scale %s -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config single_pointing.ms_aterm_norotate_config.txt -size %d %d %s/single_pointing.ms'%(inputs['wsclean_exec'],inputs['output_path'],inputs['weight'],inputs['cell'],int(inputs['size']),int(inputs['size']),inputs['output_path']))

	if inputs['mosaic'] == "False":
		commands.append('%s %s/single_pointing-image-pb.fits'%(inputs['rms_exec'],inputs['output_path']))
	else:
		commands.append('%s simulations/fit_pb.py'%(inputs['CASA_exec']))
		commands.append('%s simulations/generate_mosaic_pointings.py'%(inputs['CASA_exec']))
		commands.append('%s simulations/make_measurement_set.py mosaic simulator_inputs.txt'%(inputs['stimela_exec']))

	with open('job_%s.%s'%(step,params['job_manager']), 'a') as filehandle:
		for listitem in commands:
			filehandle.write('%s\n' % listitem)

if part==2:
	if inputs['mosaic'] == "True":
		commands = []
		step = 'mosaic'
		write_hpc_headers(step,params)

		commands.append('array=(%s/mosaic_*.ms)'%inputs['output_path'])
		commands.append('len=${#array[@]}')
		commands.append('a=$SLURM_ARRAY_TASK_ID')
		## Add noise to all ms
		commands.append('%s simulations/add_noise_hetero.py ${array[$a]} %d %.3f %s %s'%(inputs['CASA_exec'],int(inputs['size']),float(inputs['time_multiplier']),band,inputs['cell']))

		## Make all a terms
		commands.append('%s simulations/generate_pb_aterms.py ${array[$a]} 0 0 0 %s'%(inputs['CASA_exec'],band))

		## Unzip a terms
		commands.append('gunzip -f ${array[$a]}\"_pb_flat_norotate.fits.gz\"')

		## Make images
		commands.append('%s -name %s/${array[$a]}_IM -no-update-model-required --aterm-kernel-size 157 -weight %s -scale %s -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config ${array[$a]}_aterm_norotate_config.txt -size %d %d ${array[$a]}'%(inputs['wsclean_exec'],inputs['output_path'],inputs['weight'],inputs['cell'],int(inputs['size']),int(inputs['size'])))

		## Convert to casa ims
		commands.append('%s simulations/convert_fits_to_casa.py ${array[$a]}'%inputs['CASA_exec'])

		with open('job_%s.%s'%(step,params['job_manager']), 'a') as filehandle:
			for listitem in commands:
				filehandle.write('%s\n' % listitem)
		
		commands = []
		step = 'make_image'
		write_hpc_headers(step,params)

		## Make mosaic
		commands.append('%s simulations/make_mosaic.py'%inputs['CASA_exec'])

		## Make rms map
		commands.append('%s %s/mosaic.linmos.fits'%(inputs['rms_exec'],inputs['output_path']))

		with open('job_%s.%s'%(step,params['job_manager']), 'a') as filehandle:
			for listitem in commands:
				filehandle.write('%s\n' % listitem)
