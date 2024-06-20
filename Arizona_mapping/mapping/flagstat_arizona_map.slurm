#!/bin/bash

#SBATCH --job-name flagstat
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e bwa_errs_outs/err_flagstat_%A_%a.err
#SBATCH -o bwa_errs_outs/std_flagstat_%A_%a.out
#SBATCH --array 1-41



# load modules necessary
module load gcc/12.2.0 samtools/1.16.1

# reference genome path
REF=/project/getpop/ref/GCA_022577455.1_rAriEle1.0.p_genomic.fna

# bam files
BAM_DIR=/project/getpop/bwa_arizona

# output director
OUT_DIR=/project/getpop/bwa_arizona/flagstat




# Set working directory to where the sorted bam files are
cd $BAM_DIR
#array of files
for x in *.bam
do   
  bams=(${bams[@]} "${x}")
done

sample=${bams[($SLURM_ARRAY_TASK_ID-1)]}
#output file name
name=$(echo $sample| cut -d '_' -f 1)

#run
samtools flagstat -@ 12 $sample > $OUT_DIR/${name}_flagstat.out





