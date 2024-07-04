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
#SBATCH -e /project/getpop/scripts/errs_outs/bwa_errs_outs/err_flagstat_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/bwa_errs_outs/std_flagstat_%A_%a.out
#SBATCH --array 1-41



# load modules necessary
module load gcc/12.2.0 samtools/1.16.1

REF=/project/getpop/ref/ratsnake_genome_files/ratsnake_scaffold/P_alleghaniensis_concatenate_ragtag.scaffolds.fasta

OUT_DIR=/project/getpop/bwa_ratsnake/flagstat

mkdir -p $OUT_DIR

# Set working directory to where the sorted reads are

cd /project/getpop/bwa_ratsnake
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





