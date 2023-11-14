#!/bin/bash

#SBATCH --job-name Trim2
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24   
#SBATCH --mem=36G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=seanharrington256@gmail.com
#SBATCH -e err_trim2_%A_%a.err
#SBATCH -o std_trim2_%A_%a.out
#SBATCH --array=1-41

#make sure cpus match threds

# Try out trimming using adapters.fa from bbduk

#Load modules
module load swset/2018.05 gcc/7.3.0 trimmomatic/0.36

#Set working directory
cd /project/getpop/raw_data

#Array
for x in *
do
	dirs=(${dirs[@]} "${x}")
done

#Call directory, -1 will let you use array number due to count starting at 0
cd ${dirs[($SLURM_ARRAY_TASK_ID-1)]}

#Lable reads 
FORWARD=*_R1_*
REVERSE=*_R2_*
ADAPTERS='/project/getpop/adapters.fa'

#make file location
SAMPLE=$(echo ${dirs[($SLURM_ARRAY_TASK_ID-1)]} | cut -f 1 -d '_')
mkdir /project/getpop/trim_out_bbduk_adapters/trim_read_$SAMPLE

OUT_DIR=/project/getpop/trim_out_bbduk_adapters/trim_read_$SAMPLE

o_for="$(basename $FORWARD)"
o_rev="$(basename $REVERSE)"

#script for trimmed reads
trimmomatic PE -threads 24 $FORWARD $REVERSE\
    $OUT_DIR/paired-$o_for $OUT_DIR/unpaired-$o_for\
    $OUT_DIR/paired-$o_rev $OUT_DIR/unpaired-$o_rev\
    ILLUMINACLIP:$ADAPTERS:2:30:10:8:true\
    LEADING:10 TRAILING:10\
    SLIDINGWINDOW:4:20 MINLEN:36
