#!/bin/bash

#SBATCH --job-name fastqc_rawdata
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hannahvantassel44@gmail.com
#SBATCH -e err_fastqc_%A_%a.err
#SBATCH -o std_fastqc_%A_%a.out
#SBATCH --array=1-46

# load up modules
module load swset/2018.05 gcc/7.3.0 fastqc/0.11.7

# set working directory
cd /project/getpop/raw_data/

#Make array containing all the directories with data
for x in *
do
	dirs=(${dirs[@]} "${x}")
done

#Call directory, -1 will let you use array number due to count starting at 0
cd ${dirs[($SLURM_ARRAY_TASK_ID-1)]}


# run fastqc on all files, output to different directory
fastqc *_R1_* *_R2_* -t 2 -o /project/getpop/fastqc_out

