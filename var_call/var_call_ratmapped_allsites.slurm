#!/bin/bash

#SBATCH --job-name ratmap_var_call
#SBATCH -A inbreh
#SBATCH -t 0-72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_out_vcf/err_varcallA_ratmap_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/err_out_vcf/std_varcallA_ratmap_%A_%a.out
#SBATCH --array=1-199



# Variant calling all samples at once


#load modules
module load gcc/12.2.0 bcftools/1.16

## Create directory path 
sorted=/project/getpop/bwa_ratsnake/rmd

#reference
REF=/project/getpop/ref/ratsnake_genome_files/ratsnake_scaffold/P_alleghaniensis_concatenate_ragtag.scaffolds.fasta

# Get the region of interest from the bed file
scaf=`cat /project/getpop/ref/ratsnake_genome_files/ratsnake_scaffold/P_alleghaniensis_concatenate_ragtag.scaffolds.bed | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`

echo "started scaffold $scaf"

## Create some output folders
OUTDIR=/project/getpop/vcf_rat_map_all
mkdir -p $OUTDIR

# Output path
OUT=${OUTDIR}/ratmap_allsites

# Go to target directory
cd $sorted

#run mpileup
bcftools mpileup --threads 12 -A -a DP,AD --regions $scaf -f $REF *.rmd.bam | bcftools call -f GQ,GP -m --threads 12 -O z -o ${OUT}_$scaf.vcf.gz

echo "completed scaffold $scaf"