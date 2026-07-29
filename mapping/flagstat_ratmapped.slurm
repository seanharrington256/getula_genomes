#!/bin/bash

#SBATCH --job-name flagstat
#SBATCH -A inbreh
#SBATCH -t 0-36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_flagstat_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_flagstat_%A_%a.out
#SBATCH --array 1-8



# load modules necessary
module load gcc/14.2.0 samtools/1.20

REF=/project/getpop/ref/ratsnake_genome_files/ratsnake_scaffold/P_alleghaniensis_concatenate_ragtag.scaffolds.fasta

OUT_DIR=/project/getpop/bwa_ratsnake/flagstat_admera

mkdir -p $OUT_DIR

# Set working directory to where the mapped reads are
cd /project/getpop/admera_bwa_ratsnake

#array of files
bams=(*.bam)


sample=${bams[($SLURM_ARRAY_TASK_ID-1)]}
#output file name
name="${sample%_sort.bam}"



#run
samtools flagstat -@ 12 $sample > $OUT_DIR/${name}_flagstat.out





