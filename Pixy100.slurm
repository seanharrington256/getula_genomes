#!/bin/bash

#SBATCH --job-name Pixy
#SBATCH -A getpop
#SBATCH -t 0-2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dahuja1@uwyo.edu
#SBATCH -e errs_outs/err_Pixy_%A_%a.err
#SBATCH -o errs_outs/std_Pixy_%A_%a.out
#SBATCH --array=1-18


# load modules necessary
module load miniconda3/23.11.0

# Activate pixy environment
conda activate pixy

# Set working directory to where the index files and trimmed_reads directory are
cd /project/getpop/vcf_allsites/filt_scafs_vcfs

OUT_DIR=/project/getpop/pixy_out100

# define population file
popfile=/project/getpop/scripts/pixy_popfiles/INDS_POP_EAST_WEST.txt


# use a loop to find all the files with pattern filtered_sort_NC_*.vcf.gz
#     and assign them to a bash array

for x in filtered_sort_NC_*.vcf.gz
do   
  infiles=(${infiles[@]} "${x}")
done


## For whichever SLURM_ARRAY_TASK_ID index a job is in, get the sample 
## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero
##   I think it's less confusing to subtract 1 here than to remember to do it when 
##   specifying the number of jobs for the array

sample=${infiles[($SLURM_ARRAY_TASK_ID-1)]}

filepre="${sample%%.*}"

# run pixy
pixy --stats pi fst dxy \
--vcf $sample \
--population $popfile \
--window_size 100000 \
--n_cores 8 \
--output_folder $OUT_DIR \
--output_prefix $filepre 


