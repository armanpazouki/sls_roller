#!/bin/bash

#SBATCH --job-name=arrayJob
#SBATCH --array=1-660
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH -o /dev/null
##SBATCH --nodelist=euler[3-21]
#SBATCH --exclude=euler[23-53]

##cd $PBS_O_WORKDIR

###KEEP THE NUMBERS THE SAME WT and ppn
###cat /etc/hostname > /home/pazouki/Renders/vehicle/a3/run_host/$PBS_ARRAYID.runhost
### povray +WT1 +KFI1 +KFF4000 +SF$PBS_ARRAYID +EF$PBS_ARRAYID renderZ_messedUp.ini
povray +KFI1 +KFF4000 +SF$SLURM_ARRAY_TASK_ID +EF$SLURM_ARRAY_TASK_ID render_batch.ini
##povray +SF$SLURM_ARRAY_TASK_ID +EF$SLURM_ARRAY_TASK_ID render_batch.ini
