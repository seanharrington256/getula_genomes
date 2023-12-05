#!/bin/bash

#SBATCH --job-name vcf_comb_stats
#SBATCH -A inbreh
#SBATCH -p teton
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e errs_outs/err_out_vcf/err_vcfstats_%A.err
#SBATCH -o errs_outs/err_out_vcf/std_vcfstats_%A.out

#load packages
module load gcc/12.2.0 bcftools/1.16 vcftools/0.1.16 r/4.2.2


# combines the pre-filtered and post-filtered vcfs from NC scaffolds
#      gets some statistics about the vcf files
#      then runs `r_vcf_stats_cmdline.R` to make plots.


# number of unfiltered variants - removing this for now
# bcftools view -H full_only100_thin20k.vcf.gz | wc -l #may take some time


# set some paths
PREFILT="/project/getpop/vcf_allsites/sort_NC_*.vcf.gz"
PREFILTCOMB="/project/getpop/vcf_allsites/sort_NC_allsites_combined.vcf.gz"

POSTFILT="/project/getpop/vcf_allsites/filt_scafs_vcfs/filtered_sort_NC*.vcf.gz"
POSTFILTCOMB="/project/getpop/vcf_allsites/filt_scafs_vcfs/filt_sort_NC_allsites_combined.vcf.gz"

VCFSTATSDIR="/project/getpop/vcf_stats/"

################# COMBINE FILES

# COMBINE PRE-FILTERED
bcftools concat -Oz -o $PREFILTCOMB $PREFILT

# COMBINE POST-FILTERED
bcftools concat -Oz -o $POSTFILTCOMB $POSTFILT



############# GET STATS

# Loop to create vcf stats for each combined file



for VCFFile in "$PREFILTCOMB" "$POSTFILTCOMB"; do

	VCF_W_EXT=$(basename "$VCFFile")

	VCFNAME="${VCF_W_EXT%.vcf.gz}"
	
	OUT=${VCFSTATSDIR}${VCFNAME}
	mkdir -p $OUT
	OUTFILE=vcf_stats_${VCFNAME}
	
	# index vcf file
	bcftools index $VCFFile
	
	
	##calc allele frequency
	vcftools --gzvcf $VCFFile --freq2 --max-alleles 2 --out $OUT/$OUTFILE
	##calc mean depth/individual
	vcftools --gzvcf $VCFFile --depth --out $OUT/$OUTFILE
	
	##calc mean depth/site
	vcftools --gzvcf $VCFFile --site-mean-depth --out $OUT/$OUTFILE
	
	##site quality
	vcftools --gzvcf $VCFFile --site-quality --out $OUT/$OUTFILE
	##proportion of missing/individual
	vcftools --gzvcf $VCFFile --missing-indv --out $OUT/$OUTFILE
	
	##proportion of missing/site
	vcftools --gzvcf $VCFFile --missing-site --out $OUT/$OUTFILE
	
	##calc heterozygosity & inbreeding coef/individual
	vcftools --gzvcf $VCFFile --het --out $OUT/$OUTFILE
	
	
	# Use R script to make plots:
	Rscript r_vcf_stats_cmdline.R $OUT $OUTFILE
done


