#!/bin/bash

#SBATCH --job-name rmd
#SBATCH -A inbreh
#SBATCH -t 0-72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_rmd_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_rmd_%A_%a.out
#SBATCH --array 1-8


# Load modules
module load gcc/14.2.0 picard/3.2.0 samtools/1.20

# Set up output directory variable
OUTDIR=/project/getpop/bwa_ratsnake/rmd_admera

mkdir -p $OUTDIR

# Set working dir to where the bam files are:
cd /project/getpop/admera_bwa_ratsnake

#array of files
bams=(*.bam)


# Get the file name of each sample
sample=${bams[($SLURM_ARRAY_TASK_ID-1)]}
# get individual name only from that file name
name=${sample%_sort.bam}


#add readgroups
picard AddOrReplaceReadGroups \
    -I $sample \
    -O $OUTDIR/${name}.rg.bam \
    -RGID id \
    -RGLB lib1 \
    -RGPL illumina \
    -RGPU unit1 \
    -RGSM $name



#run picard MarkDuplicates
picard \
 MarkDuplicates -REMOVE_DUPLICATES true \
 -ASSUME_SORTED true -VALIDATION_STRINGENCY SILENT \
 -INPUT $OUTDIR/${name}.rg.bam \
 -OUTPUT $OUTDIR/${name}.rmd.bam \
 -METRICS_FILE $OUTDIR/${name}.rmd.bam.metrics
 
 
samtools index $OUTDIR/${name}.rmd.bam
