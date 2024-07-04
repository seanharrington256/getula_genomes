#!/bin/bash

#SBATCH --job-name pixy
#SBATCH -A getpop
#SBATCH -t 0-4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/getpop/scripts/errs_outs/err_Pixy_%A_%a.err
#SBATCH -o /project/getpop/scripts/errs_outs/std_Pixy_%A_%a.out

# load modules necessary
module load miniconda3/23.11.0

# Activate pixy environment
conda activate pixy

# Set working directory to where the index files and trimmed_reads directory are
cd /project/getpop/vcf_allsites_Arizona_map

OUT_DIR=/project/getpop/pixy_out_AriMap/pixy_out10_combined

mkdir -p $OUT_DIR

# define population file
popfile=/project/getpop/scripts/gen_scan/pixy_popfiles/INDS_POP_EAST_WEST.txt

sample=filtered_getula_arimap_AllSites_full_genome.vcf.gz

# filepre="${sample%%.*}"

filepre="arimap10k"

tabix $sample

# run pixy
pixy --stats pi fst dxy \
--vcf $sample \
--population $popfile \
--window_size 10000 \
--n_cores 16 \
--output_folder $OUT_DIR \
--output_prefix $filepre 


