#!/bin/bash

#SBATCH --job-name var_call
#SBATCH -A inbreh
#SBATCH -t 0-72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_varcallArizAllSite_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_varcallArizAllSite_%A_%a.out
#SBATCH --array=1-140



# Variant calling all samples at once


#load modules
module load gcc/12.2.0 bcftools/1.16

## Create directory path to sorted bams
sorted=/project/getpop/bwa_arizona/

#reference
REF=/project/getpop/ref/GCA_022577455.1_rAriEle1.0.p_genomic.fna


# Get the region of interest from the bed file
scaf=`cat /project/getpop/ref/Arizona_ref_bed.bed | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`


## Create some output folders
OUTDIR=/project/getpop/vcf_allsites_Arizona_map_NORMD/
mkdir -p $OUTDIR



echo "started scaffold $scaf"


# Go to target directory
cd $sorted

#run mpileup
bcftools mpileup --threads 12 -A -a DP,AD --regions $scaf -f $REF *.bam | bcftools call -f GQ,GP -m -Oz -f GQ --threads 12 -o ${OUTDIR}${scaf}_allsites_unfiltNORMD.vcf.gz

echo "completed scaffold $scaf"
