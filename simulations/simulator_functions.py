import re, os, inspect, sys
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from astropy import wcs
import astropy.units as u
from astropy.io import fits
from matplotlib.colors import SymLogNorm
import pandas as pd


def write_hpc_headers(step,params):
	func_name = inspect.stack()[0][3]
	print(params)
	hpc_opts = {}
	hpc_opts['job_manager'] = params['job_manager']
	hpc_opts['job_name'] = 'vlbisched_%s'%step
	hpc_opts['email_progress'] = params["email_progress"] 
	hpc_opts['hpc_account'] = params['HPC_project_code']
	hpc_opts['error'] = step

	if ((hpc_opts['job_manager'] == 'pbs')|(hpc_opts['job_manager'] == 'bash')|(hpc_opts['job_manager'] == 'slurm')):
		pass
	#else:
		#sys.exit()

	for i in ['partition','walltime','nodetype','nodes','cpus','mpiprocs']:
		hpc_opts[i] = params['%s'%i]


	hpc_dict = {'slurm':{
					 'partition'     :'#SBATCH --partition=%s'%hpc_opts['partition'],
					 'nodetype'      :'',
					 'cpus'          :'#SBATCH --tasks-per-node %s'%hpc_opts['cpus'], 
					 'nodes'         :'#SBATCH -N %s-%s'%(hpc_opts['nodes'],hpc_opts['nodes']),
					 'mpiprocs'      :'', 
					 'walltime'      :'#SBATCH --time=%s'%hpc_opts['walltime'],
					 'job_name'      :'#SBATCH -J %s'%hpc_opts['job_name'],
					 'hpc_account'   :'#SBATCH --account %s'%hpc_opts['hpc_account'],
					 'email_progress':'#SBATCH --mail-type=BEGIN,END,FAIL\n#SBATCH --mail-user=%s'%hpc_opts['email_progress'],
					 'error':'#SBATCH -o logs/%s.sh.stdout.log\n#SBATCH -e logs/%s.sh.stderr.log'%(hpc_opts['error'],hpc_opts['error'])
					},
				'pbs':{
					 'partition'     :'#PBS -q %s'%hpc_opts['partition'],
					 'nodetype'      :'',
					 'cpus'          :'#PBS -l select=%s:ncpus=%s:mpiprocs=%s:nodetype=%s'%(hpc_opts['nodes'],hpc_opts['cpus'],hpc_opts['mpiprocs'],hpc_opts['nodetype']), 
					 'nodes'         :'',
					 'mpiprocs'      :'', 
					 'walltime'      :'#PBS -l walltime=%s'%hpc_opts['walltime'],
					 'job_name'      :'#PBS -N %s'%hpc_opts['job_name'],
					 'hpc_account'   :'#PBS -P %s'%hpc_opts['hpc_account'],
					 'email_progress':'#PBS -m abe -M %s'%hpc_opts['email_progress'],
					 'error':'#PBS -o logs/%s.sh.stdout.log\n#PBS -e logs/%s.sh.stderr.log'%(hpc_opts['error'],hpc_opts['error'])
					},
				'bash':{
					 'partition'     :'',
					 'nodetype'      :'',
					 'cpus'          :'', 
					 'nodes'         :'',
					 'mpiprocs'      :'', 
					 'walltime'      :'',
					 'job_name'      :'',
					 'hpc_account'   :'',
					 'email_progress':'',
					 'error':''
					}
				}

	hpc_header= ['#!/bin/bash']

	if step == 'run_mosaic_sims':
		file = open("%s/target_files.txt"%params['global']['cwd'], "r")
		nonempty_lines = [line.strip("\n") for line in file if line != "\n"]
		line_count = len(nonempty_lines)
		file.close()
		if params[step]['hpc_options']['max_jobs'] == -1:
			tasks = '0-'+str(line_count-1)
		else:
			if (line_count-1) > params[step]['hpc_options']['max_jobs']:
				tasks = '0-'+str(line_count-1)+'%'+str(params[step]['hpc_options']['max_jobs'])
			else:
				tasks = '0-'+str(line_count-1)
		hpc_dict['slurm']['array_job'] = '#SBATCH --array='+tasks
		hpc_dict['pbs']['array_job'] = '#PBS -t '+tasks
		hpc_dict['bash']['array_job'] = ''
		hpc_opts['array_job'] = -1

	hpc_job = hpc_opts['job_manager']
	for i in hpc_opts.keys():
		if i != 'job_manager':
			if hpc_opts[i] != '':
				if hpc_dict[hpc_opts['job_manager']][i] !='':
					hpc_header.append(hpc_dict[hpc_job][i])


	with open('job_%s.%s'%(step,hpc_job), 'w') as filehandle:
		for listitem in hpc_header:
			filehandle.write('%s\n' % listitem)