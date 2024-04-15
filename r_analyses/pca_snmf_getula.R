# R script to run PCA and snmf from the LEA package on getula genomes

# load up stuff on cluster:
# module load gcc r/4.2.2 rstudio
# module load udunits/2.2.28 gdal/3.6.1

library(LEA)
library(ggplot2)
library(scatterpie)
library(sf)
library(maps)
library(mapdata)
library(rworldmap)



# set up some paths
lfmm_dir <- "/pfs/tc1/project/getpop/lfmm/"
# ped_file <- "filtered_ratmapAll_NConly.ped"
coords_file <- "/project/getpop/metadata/getula_genome_coords.csv"

# ## test on a smaller chromosome at first:
ped_file <- "rat_map_all_NC_045557.ped"


setwd(lfmm_dir)
ind_names <- readLines("ind_names_lfmm.txt")

#### Some overall setup for mapping and plotting

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")

## Read in coordinates for plotting farther down
coords<-read.csv(coords_file, header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't


################################################# map data
sf_use_s2(FALSE) # turn off sphere geometry in sf - doesn't play nicely with the shape files I'm using

# Set up boundaries for map plotting:
xmin <- -125
xmax <- -65
ymin <- 25
ymax <- 50


# Define Albers Equal Area Conic projection
target_crs <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Set up basemap of eastern NA
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) # Get US states data
mexico <- st_as_sf(map("worldHires", "Mexico", plot = FALSE, fill = TRUE)) # Get Mexico data
canada <- st_as_sf(map("worldHires", "Canada", plot = FALSE, fill = TRUE)) # Get Canada data
to_map <- rbind(states, mexico, canada) # Combine the datasets

# crop the base map polygon - based on this: https://datascience.blog.wzb.eu/2019/04/30/zooming-in-on-maps-with-sf-and-ggplot2/
to_map_val <- st_make_valid(to_map)
cropped_map <- st_crop(to_map_val, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

# make a ggplot object with the basemap
basemap <- ggplot() +
  geom_sf(data = cropped_map, fill = "gray95", color = "black") +
  theme_minimal()








# convert the ped file
ped2geno(ped_file, force = FALSE)
ped2lfmm(ped_file, force = FALSE)
lfmm_file <- gsub("ped$", "lfmm", ped_file)
geno_file <- gsub("ped$", "geno", ped_file)

# When trying to a run a PCA below, I ran into this issue for one chromosome (NC_045557):
#   Error: SNP 8 is constant among individuals.
#   this isnt *actually* constant, but is heterozygous in all individuals
#   this presumably also carries no info

# Filter the LFMM file to remove any loci where there are only heterozygotes
lfmm_dat <- read.table(lfmm_file)

# Function to check if a column has only one unique value, not including missing data
more_one_uniq <- function(x) {
  length(unique(x[x !=9])) > 1 # is there more than 1 unique value, not including 9?
}

# Apply the function across columns and subset the dataframe
lfmm_dat_filtered <- lfmm_dat[, apply(lfmm_dat, 2, more_one_uniq)]
lfmm_filt_file <- paste0("pca_filt_", lfmm_file)
write.table(lfmm_dat_filtered, file = lfmm_filt_file, row.names = FALSE, col.names = FALSE)


## Run a PCA:
pc = pca(lfmm_filt_file, scale = TRUE)
tw = tracy.widom(pc)
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
# suggests 2-4 pops, maybe up to 5?

# PC1-PC2 plot - bad plot, can clean this up later
plot(pc$projections)

# Run snmf
# Run sNMF using 1 to 10 ancestral populations and evaluate the fit of different k values to the data using cross entropy criterion
# before running snmf, check if it's already been run and then just load if it has
if(dir.exists(gsub("geno", "snmf", basename(geno_file)))){
  obj.at<-load.snmfProject(gsub("geno", "snmfProject", basename(geno_file))) # if it has, just load up the results
}else{ # otherwise, run sNMF
  obj.at <- snmf(input.file = geno_file,  # input file is the .geno format file. We set up the path to this above
                 K = 1:10, # we will test for k=1 through 10
                 ploidy = 2, 
                 entropy = T, # use the cross entropy criterion for assessing the best k value
                 repetitions = 10, # Run 10 independent replicate analyses
                 CPU = 6, 
                 project = "new", tolerance = 0.00001, iterations = 500)
}

pdf(file = "crossentropy.pdf", width = 8, height=5)
plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
dev.off()













