#!/bin/bash

#SBATCH --job-name Pixy
#SBATCH -A getpop
#SBATCH -t 0-2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dahuja1@uwyo.edu
#SBATCH -e errs_outs/err_Pixy_%A_%a.err
#SBATCH -o errs_outs/std_Pixy_%A_%a.out


# load modules necessary
module load miniconda3/23.11.0

# Activate pixy environment
conda activate pixy

# Set working directory to where the index files and trimmed_reads directory are
cd /project/getpop/vcf_allsites/filt_scafs_vcfs

OUT_DIR=/project/getpop/pixy_out/pixy_out100_combined

# define population file
popfile=/project/getpop/scripts/pixy_popfiles/INDS_POP_EAST_WEST.txt

sample=filt_sort_NC_allsites_combined.vcf.gz

filepre="${sample%%.*}"

# run pixy
pixy --stats pi fst dxy \
--vcf $sample \
--population $popfile \
--window_size 100000 \
--n_cores 16 \
--output_folder $OUT_DIR \
--output_prefix $filepre 


