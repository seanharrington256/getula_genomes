#!/bin/bash

#SBATCH --job-name vcfstats
#SBATCH -A inbreh
#SBATCH -t 0-8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e err_out_vcf/err_filthand_%A.err
#SBATCH -o err_out_vcf/std_filthand_%A.out

#load packages
module load gcc/12.2.0 bcftools/1.16 vcftools/0.1.16

#move to directory 
cd /project/getpop/vcf_rat_map_all




# number of unfiltered variants - removing this for now
# bcftools view -H full_only100_thin20k.vcf.gz | wc -l #may take some time


# Set up files

VCFFile=/project/getpop/vcf_rat_map_all/filtered_ratmapAll_full_genome.vcf.gz
OUT=/project/getpop/vcf_stats/vcf_stats_ratmapAllcallFilt
mkdir -p $OUT
OUTFILE=vcf_stats_ratmapAllcallFilt

# index vcf file
bcftools index $VCFFile



##calc allele frequency
echo "allele frequency"
vcftools --gzvcf $VCFFile --freq2 --max-alleles 2 --out $OUT/$OUTFILE
##calc mean depth/individual
echo "mean depth/individual"
vcftools --gzvcf $VCFFile --depth --out $OUT/$OUTFILE

##calc mean depth/site
echo "mean depth/site"
vcftools --gzvcf $VCFFile --site-mean-depth --out $OUT/$OUTFILE

##site quality
echo "site qual"
vcftools --gzvcf $VCFFile --site-quality --out $OUT/$OUTFILE
##proportion of missing/individual
echo "proportion of missing/individual"
vcftools --gzvcf $VCFFile --missing-indv --out $OUT/$OUTFILE

##proportion of missing/site
echo "proportion of missing/site"
vcftools --gzvcf $VCFFile --missing-site --out $OUT/$OUTFILE

##calc heterozygosity & inbreeding coef/individual
echo "heterozygosity & inbreeding coef/individual"
vcftools --gzvcf $VCFFile --het --out $OUT/$OUTFILE
