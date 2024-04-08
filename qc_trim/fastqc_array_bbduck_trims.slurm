#!/bin/bash

#SBATCH --job-name fastqc_trimmed_bbduck
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-16:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannahvantassel44@gmail.com
#SBATCH -e err_trimbbduck_fastqc_%A_%a.err
#SBATCH -o std_trimbbduck_fastqc_%A_%a.out
#SBATCH --array=1-41

# load up modules
module load swset/2018.05 gcc/7.3.0 fastqc/0.11.7

# set working directory
cd /project/getpop/trim_out_bbduk_adapters

#Make array containing all the directories with data
for x in *
do
	dirs=(${dirs[@]} "${x}")
done

#Call directory, -1 will let you use array number due to count starting at 0
cd ${dirs[($SLURM_ARRAY_TASK_ID-1)]}


# run fastqc on all files, output to different directory
fastqc paired*_R1_* paired*_R2_* unpaired*_R1_* unpaired*_R2_* -t 4 -o /project/getpop/fastqc_out_trim_bbduck/

