#!/bin/bash

#SBATCH --job-name fastqc_trimmed
#SBATCH -A getpop
#SBATCH -t 0-2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_trimfastqc_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_trimfastqc_%A_%a.out
#SBATCH --array=1-8

# load up modules
module load fastqc/0.12.1

# set working directory
cd /project/getpop/trim_out_admera


# list out the directories
dirs=(*/)


#Call directory, -1 will let you use array number due to count starting at 0
cd ${dirs[($SLURM_ARRAY_TASK_ID-1)]}


# run fastqc on all files, output to different directory
fastqc paired*_R1_* paired*_R2_* -t 2 -o /project/getpop/fastqc/fastqc_out_admera_trimmed

