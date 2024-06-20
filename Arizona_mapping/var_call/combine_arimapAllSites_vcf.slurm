#!/bin/bash

#SBATCH --job-name combinevcfs
#SBATCH -A inbreh
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_combinevcfs_%A.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_combinevcfs_%A.out


#load packages
module load gcc/12.2.0 bcftools/1.16

#move to directory 
cd /project/getpop/vcf_allsites_Arizona_map


# combine
bcftools concat -Oz -o getula_arimap_AllSites_full_genome.vcf.gz sort_*.vcf.gz

