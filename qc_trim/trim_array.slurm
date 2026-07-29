#!/bin/bash

#SBATCH --job-name trim
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24   
#SBATCH --mem=36G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_trim2_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_trim2_%A_%a.out
#SBATCH --array=1-8


# Try out trimming using adapters.fa from bbduk

#Load modules
module load trimmomatic/0.39

#Set working directory
cd /project/getpop/admera_raw_data

# Array of directories
dirs=(*/)


SAMPLE=$(echo ${dirs[($SLURM_ARRAY_TASK_ID-1)]})

#Call directory, -1 will let you use array number due to count starting at 0
cd $SAMPLE

#Lable reads 
FORWARD=*_R1_*
REVERSE=*_R2_*
ADAPTERS='/project/getpop/adapters.fa'

#make output directory location
OUT_DIR=/project/getpop/trim_out_admera/trim_read_$SAMPLE
mkdir -p $OUT_DIR

o_for="$(basename $FORWARD)"
o_rev="$(basename $REVERSE)"

#script for trimmed reads
trimmomatic PE -threads $SLURM_CPUS_PER_TASK $FORWARD $REVERSE\
    $OUT_DIR/paired-$o_for $OUT_DIR/unpaired-$o_for\
    $OUT_DIR/paired-$o_rev $OUT_DIR/unpaired-$o_rev\
    ILLUMINACLIP:$ADAPTERS:2:30:10:8:true\
    LEADING:10 TRAILING:10\
    SLIDINGWINDOW:4:20 MINLEN:36
