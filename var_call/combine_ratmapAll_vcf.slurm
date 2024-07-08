#!/bin/bash

#SBATCH --job-name combinevcfs
#SBATCH -A inbreh
#SBATCH -t 0-20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_out_vcf/err_combinevcfs_%A.err
#SBATCH -o /project/getpop/scripts/errs_outs/err_out_vcf/std_combinevcfs_%A.out


#load packages
module load gcc/12.2.0 bcftools/1.16

#move to directory 
cd /project/getpop/vcf_rat_map_all


# combine
bcftools concat -Oz -o ratmapAllsites_genome.vcf.gz sort_ratmap_allsites*.vcf.gz

bcftools concat -Oz -o NC_scafsAllsites_genome.vcf.gz sort_ratmap_allsites_NC_*.vcf.gz