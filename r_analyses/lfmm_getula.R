## Script to run LFMM on getula genomes and identify any loci that might
##    be under environmental selection

# load necessary libraries
library(LEA)



# set up some paths
lfmm_dir <- "/pfs/tc1/project/getpop/lfmm/"
ped_file <- "filtered_ratmapAll_NConly.ped"
env_file <- "/pfs/tc1/project/getpop/metadata/enviro_data/getula_coords_red_enviro_vars.csv"

# ## test on a smaller chromosome at first:
# ped_file <- "filtered_sort_rat_map_all_NC_045557.1_RagTag.ped"

setwd(lfmm_dir)
ind_names <- readLines("ind_names_lfmm.txt")


# convert the ped file
ped2lfmm(ped_file, force = FALSE)


# read in environmental data:
env_data <- read.csv(env_file)

# variables to start with: AnnMeanT_B1 (annual mean temp), AnnPrecip_B12 (annual precipitation)



