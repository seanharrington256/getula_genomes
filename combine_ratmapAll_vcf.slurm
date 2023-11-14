#!/bin/bash

#SBATCH --job-name combinevcfs
#SBATCH -A inbreh
#SBATCH -t 0-15:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e err_out_vcferr_combinevcfs_%A.err
#SBATCH -o err_out_vcfstd_combinevcfs_%A.out


#load packages
module load gcc/12.2.0 bcftools/1.16

#move to directory 
cd /project/getpop/vcf_rat_map_all


# combine
bcftools concat -Oz -o ratmapAll_full_genome.vcf.gz sort_rat_map_all*.vcf.gz

