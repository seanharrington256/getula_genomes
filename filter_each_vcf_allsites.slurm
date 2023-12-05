#!/bin/bash

#SBATCH --job-name vcf_filt
#SBATCH -A inbreh
#SBATCH -t 0-64:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e errs_outs/err_out_vcf/err_filt_%A_%a.err
#SBATCH -o errs_outs/err_out_vcf/std_filt_%A_%a.out
#SBATCH --array=1-199


#load packages
module load gcc/12.2.0 bcftools/1.16 vcftools/0.1.16

#go to directory
cd /project/getpop/vcf_allsites

# set up the inputs:
#VCF_IN=ratmapAll_full_genome.vcf.gz
#VCF_OUT=filtered_$VCF_IN

mkdir -p filt_scafs_vcfs

#Array
#make array
for x in sort*.vcf.gz
do   
  vcf_all=(${vcf_all[@]} "${x}")
done

vcf=${vcf_all[($SLURM_ARRAY_TASK_ID-1)]}

# name=$(echo $vcf | cut -d '_' -f 3 | cut -d '.' -f 1)

##Apply Filters to VCF
VCF_IN=$vcf
VCF_OUT="filt_scafs_vcfs/filtered_$vcf"

#filters
MAF=0  # set this to zero for correct calculation of pi & dxy
MISS=0.8
QUAL=30
MIN_DEPTH=10
MIN_ALLELES=0
MAX_ALLELES=2


echo "start filtering $vcf"


#filter
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--minDP $MIN_DEPTH --min-alleles $MIN_ALLELES --max-alleles $MAX_ALLELES --recode --stdout | bgzip > \
$VCF_OUT

# index each
tabix $VCF_OUT

echo "completed filtering $vcf"

#cd $filtered_multi_vcf
#view
#bcftools view -H filtered_$name.vcf.gz | wc -l
