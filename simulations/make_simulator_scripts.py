import inspect, os, sys, json
import copy
## Python 2 will need to adjust for casa 6
import collections

filename = inspect.getframeinfo(inspect.currentframe()).filename
sys.path.append(os.path.dirname(os.path.realpath(filename)))

from simulator_functions import *

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
params['HPC_project_code'] = str(params['HPC_project_code'])
params['partition'] = str(params['partition'])
params['walltime'] = str(params['walltime'])
params['nodetype'] = str(params['nodetype'])
params['nodes'] = int(params['nodes'])
params['cpus'] = int(params['cpus'])
params['mpiprocs'] = int(params['mpiprocs'])

write_hpc_headers('test',params)

#SBATCH -N 1-1
#SBATCH --tasks-per-node 12
#SBATCH -J simmosaic
#SBATCH -m cyclic
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jack.f.radcliffe@gmail.com
#SBATCH -o logs/simmosaic.sh.stdout.log
#SBATCH -e logs/simmosaic.sh.stderr.log
#SBATCH --partition=GPU
#SBATCH --mem=50G
#SBATCH -t 01:00:00