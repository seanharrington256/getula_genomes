#!/bin/bash

#SBATCH --job-name sort_vcfs
#SBATCH -A inbreh
#SBATCH -t 0-48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_sort_ariMap_vcfs_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_sort_ariMap_vcfs_%A_%a.out
#SBATCH --array=1-140

#load packages
module load gcc/12.2.0 bcftools/1.16

cd /project/getpop/vcf_allsites_Arizona_map # set directory

#Array
for x in *allsites_unfilt.vcf.gz
do
	files=(${files[@]} "${x}")
done

vcf=${files[($SLURM_ARRAY_TASK_ID-1)]} # define the vcf file to work on

echo "started file $vcf"

# run sort
bcftools sort -Oz $vcf -o sort_$vcf

echo "completed file $vcf"

# delete the original vcf, we don't need it anymore
# rm $vcf

