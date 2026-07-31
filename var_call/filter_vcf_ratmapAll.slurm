#!/bin/bash

#SBATCH --job-name filt_get_vcf
#SBATCH -A inbreh
#SBATCH -t 0-71:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_out_vcf/err_filt_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/err_out_vcf/std_filt_%A_%a.out


#load packages
module load gcc/14.2.0 bcftools/1.20 vcftools/0.1.17

#go to directory
cd /project/getpop/vcf_rat_map_all

# set up the inputs:
VCF_IN=NC_scafsAllsites_genome.vcf.gz
VCF_OUT=filtered_$VCF_IN


# #Array   # just running this on a single VCF presently
# #make array
#   vcf_all=(sort_*.vcf.gz)

# vcf=${vcf[($SLURM_ARRAY_TASK_ID-1)]}

# name=$(echo $vcf | cut -d '_' -f 3 | cut -d '.' -f 1)

##Apply Filters to VCF
# VCF_IN=/project/getpop/sorted_multi_vcf/$vcf
# VCF_OUT=/project/getpop/filtered_multi_vcf/filtered_$name.vcf.gz

#filters
MAF=0  # set this to zero for correct calculation of pi & dxy
MISS=0.8
QUAL=30
MIN_DEPTH=7
MIN_ALLELES=2
MAX_ALLELES=2


#filter
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--minDP $MIN_DEPTH --min-alleles $MIN_ALLELES --max-alleles $MAX_ALLELES --recode --stdout | gzip -c > \
$VCF_OUT

#cd $filtered_multi_vcf
#view
#bcftools view -H filtered_$name.vcf.gz | wc -l
