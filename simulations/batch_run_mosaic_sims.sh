#!/bin/bash
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

singularity exec /idia/software/containers/casa-6.3.simg python make_itrf.py evn emerlin
singularity exec ../STIMELA_SINGULARITY_IMAGES/stimela_simms_1.2.0.img python make_measurement_set.py single evn emerlin
singularity exec /idia/software/containers/casa-6.3.simg python add_noise_hetero.py evn_emerlin.ms 2048
singularity exec /idia/software/containers/casa-6.3.simg python generate_pb_aterms.py evn_emerlin.ms 0 0 0 
gunzip "evn_emerlin.ms_pb_flat_norotate.fits.gz"
singularity exec /idia/software/containers/wsclean-gpu.simg wsclean -name evn_emerlin_nat -no-update-model-required --aterm-kernel-size 157 -weight natural -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config evn_emerlin.ms_aterm_norotate_config.txt -size 2048 2048 evn_emerlin.ms
singularity exec /idia/software/containers/wsclean-gpu.simg wsclean -name evn_emerlin_uni -no-update-model-required --aterm-kernel-size 157 -weight uniform -scale 1asec -niter 1 -mgain 0.9 -auto-threshold 0.5 -auto-mask 4 -use-idg -idg-mode hybrid -aterm-config evn_emerlin.ms_aterm_norotate_config.txt -size 2048 2048 evn_emerlin.ms