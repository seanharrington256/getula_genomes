#!/bin/bash

#SBATCH --job-name bam_idx
#SBATCH -A getpop
#SBATCH -p teton
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e bwa_errs_outs/err_bam_idx_%A_%a.err
#SBATCH -o bwa_errs_outs/std_bam_idx_%A_%a.out
#SBATCH --array=1-41


# load modules necessary
module load gcc/12.2.0 samtools/1.16.1

# Set working directory to where the index files and trimmed_reads directory are
cd /home/sharrin2/getpop/bwa_ratsnake


# use a loop to find all the files in trimmed_reads that end
#     in bam, and assign them to a bash array

for x in *bam
do   
  bams=(${bams[@]} "${x}")
done


## For whichever SLURM_ARRAY_TASK_ID index a job is in, get the file 
##     name to simplify the call below
## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero
##   I think it's less confusing to subtract 1 here than to remember to do it when 
##   specifying the number of jobs for the array

sample=${bams[($SLURM_ARRAY_TASK_ID-1)]}

# index with samtools:
samtools index -@ 12 $sample
