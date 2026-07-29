#!/bin/bash

#SBATCH --job-name BWA_rat
#SBATCH -A inbreh
#SBATCH -t 0-36:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_BWA_rat_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_BWA_rat_%A_%a.out
#SBATCH --array=1-8


# load modules necessary
module load gcc/14.2.0 bwa/0.7.17 samtools/1.20

# Set working directory to where the index files and trimmed_reads directory are
cd /project/getpop/trim_out_admera

#New ref for population loci assesment and later annotation
REF=/project/getpop/ref/ratsnake_genome_files/ratsnake_scaffold/P_alleghaniensis_concatenate_ragtag.scaffolds.fasta


OUT_DIR=/project/getpop/admera_bwa_ratsnake
mkdir -p $OUT_DIR

# use a loop to find all the files in trimmed_reads that end
#     in fastq.gz and assign them to a bash array

dirs=(*/)


## For whichever SLURM_ARRAY_TASK_ID index a job is in, get the sample 
##     name to simplify the trimmomatic call below
## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero
##   I think it's less confusing to subtract 1 here than to remember to do it when 
##   specifying the number of jobs for the array

sample=${dirs[($SLURM_ARRAY_TASK_ID-1)]}

#call directory
cd $sample

# To get the output filename, we'll do a bit of extra manipulation:
# First, strip off the stuff at the beginnig we don't want

SAMP_ID="${sample#trim_read_}"
SAMP_ID="${SAMP_ID%/}"

# strip off the .fastq.gz and add on a .sam as the output
# outname=${SAMP_ID}.sam

# Run bwa on each file using the mem algorithm, then pipe output to samtools to convert to bam and sort the bam file
bwa mem -M -t 24 $REF paired*_R1_* paired*_R2_* | samtools view - -b | samtools sort - -o $OUT_DIR/${SAMP_ID}_sort.bam

# Index the bam file
samtools index -@ 24 $OUT_DIR/${SAMP_ID}_sort.bam

