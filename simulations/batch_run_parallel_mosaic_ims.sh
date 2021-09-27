#!/bin/bash
#SBATCH -N 1-1
#SBATCH --tasks-per-node 12
#SBATCH -J simmosaic
#SBATCH -m cyclic
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jack.f.radcliffe@gmail.com
#SBATCH -o logs/simmosaic_ims.sh.stdout.log
#SBATCH -e logs/simmosaic_ims.sh.stderr.log
#SBATCH --partition=GPU
#SBATCH --array=0-170%64
#SBATCH --mem=10G
#SBATCH -t 06:00:00

array=(mosaic_ms/*.ms)
len=${#array[@]}
a=$SLURM_ARRAY_TASK_ID

singularity exec /idia/software/containers/casa-6.3.simg python add_noise_hetero.py ${array[$a]} 2048
singularity exec /idia/software/containers/casa-6.3.simg python generate_pb_aterms.py ${array[$a]} 0 0 0 
gunzip ${array[$a]}"_pb_flat_norotate.fits.gz"
singularity exec /idia/software/containers/wsclean-gpu.simg wsclean -name ${array[$a]}_image -no-update-model-required --aterm-kernel-size 157 -weight natural -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config ${array[$a]}_aterm_norotate_config.txt -size 2048 2048 ${array[$a]}
singularity exec /idia/software/containers/casa-6.3.simg python convert_fits_to_casa.py ${array[$a]}
rm casa*log