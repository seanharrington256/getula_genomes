#!/bin/bash

#SBATCH --job-name snps_vcf_filt
#SBATCH -A inbreh
#SBATCH -t 0-64:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_filt_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_filt_%A_%a.out



# Filter the SNPs down based on a minor allele count

#load packages
module load gcc/12.2.0 bcftools/1.16 vcftools/0.1.16

#go to directory
cd /project/getpop/vcf_allsites_Arizona_map

# set up the inputs:
VCF_IN=filt_getula_arimap_SNPs_MAC2_full_genome.vcf.gz
VCF_OUT=filt_getula_arimap_SNPs_MAC2_1k_THIN_full_genome.vcf.gz


#filters
MAC=2        # set this to zero for correct calculation of pi & dxy
MISS=0.8
QUAL=30
MIN_DEPTH=6
MIN_ALLELES=2    # set this to zero for correct calculation of pi & dxy
MAX_ALLELES=2    # set this to 4 for pixy (i.e., no max)


#filter
vcftools --gzvcf $VCF_IN --thin 1000 \
--remove-indels --mac $MAC --max-missing $MISS --minQ $QUAL \
--minDP $MIN_DEPTH --min-alleles $MIN_ALLELES --max-alleles $MAX_ALLELES --recode --stdout | gzip -c > \
$VCF_OUT
