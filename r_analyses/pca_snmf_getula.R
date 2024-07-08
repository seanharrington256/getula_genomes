# R script to run PCA and snmf from the LEA package on getula genomes

# load up stuff on cluster:
# module load gcc/12.2.0 r/4.4.0 rstudio
# module load udunits/2.2.28 gdal/3.6.1
# make sure base conda isn't activated/doesn't have any weird stuff in it



library(LEA)
library(ggplot2)
library(scatterpie)
library(sf)
library(maps)
library(mapdata)
library(rworldmap)



# set up some paths
lfmm_dir <- "/pfs/tc1/project/getpop/lfmm/"
ped_file <- "arimap_10k_genome.ped"
coords_file <- "/project/getpop/metadata/getula_genome_coords.csv"

# ## test on a smaller chromosome at first:


setwd(lfmm_dir)
ind_names <- readLines("ind_names_lfmm.txt")

#### Some overall setup for mapping and plotting

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "hotpink", "V4" = "purple", "V5" = "orange", "V6" = "yellow")

## Read in coordinates for plotting farther down
coords<-read.csv(coords_file, header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't
# make sure coords are in the same order as the genetic data
coords <- coords[match(ind_names, coords$number),]


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

k <- 3

ce <- cross.entropy(obj.at, K = k) 
best.run <- which.min(ce) # find the run with the lowest cross validation error

## Get the snmf Q matrix from the best run at this k
qmatrix <- Q(obj.at, K = k, run = best.run)
admix<-as.data.frame(qmatrix)

# get the coordinate and admix data into a single dataframe
for_pies <- cbind(coords, admix)

ce <- cross.entropy(obj.at, K = k) 
best.run <- which.min(ce) # find the run with the lowest cross validation error

## Get the snmf Q matrix from the best run at this k
qmatrix <- Q(obj.at, K = k, run = best.run)
admix<-as.data.frame(qmatrix)


# convert the coordinates into the target projection
for_pies_sf <- st_as_sf(for_pies, coords = c("lon", "lat"), crs = 4326) # Create an sf object with WGS84 coordinates
for_pies_sf_transformed <- st_transform(for_pies_sf, target_crs) # Transform to the target CRS
for_pies_transformed_df <- st_coordinates(for_pies_sf_transformed) %>% # Convert transformed sf object back to dataframe
  as.data.frame() %>%
  setNames(c("lon", "lat"))
for_pies_transformed_df <- cbind(for_pies[, setdiff(names(for_pies), c("lon", "lat"))], for_pies_transformed_df) # Combine with other columns from the original dataframe


# Get the right number of colors
num_cols <- length(grep("^V", colnames(for_pies_transformed_df))) # how many V columns containing proportion of cluster memberships?
colors <- colors_6[1:num_cols]

# set pie size
pie_size <- 0.7

# plot
snmf_plot <- basemap +
  geom_scatterpie(data = for_pies_transformed_df, aes(x=lon, y=lat), cols = grep("^V", colnames(for_pies_transformed_df), value = TRUE), size = 0.01, pie_scale = pie_size) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
  scale_fill_manual(values = colors) +
  guides(fill="none") + # get rid of the legend for admixture
  coord_sf(crs = target_crs) + # set target projection
  # ggtitle(sci_name) + # set title to the species
  theme(plot.title = element_text(face="bold.italic", hjust = 0.5), # format the title
        axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank())  # Remove y-axis label
snmf_plot


# Plot to pdf 
pdf(file = paste0("genome_snmf_K", k, ".pdf"), height = 13, width = 16)
print(snmf_plot)
dev.off()

# Impute the missing genotypes
impute(obj.at, geno_file, method = 'mode', K = k, run = best.run)










