## Code to extract environmental data for getula genomes


# on cluster, may need to do this to install packages: options(download.file.method="wget")
library(terra) # I need to load up my conda env and the rgdal module on the cluster to make this work
library(psych)

# set up a bunch of directories
main_dir<-"/project/getpop/metadata/enviro_data"
worldclim_dir<-"/project/getpop/metadata/enviro_data/worldclim"
bioclim_dir<-"/project/getpop/metadata/enviro_data/worldclim/wc2.1_30s_bio"
envirem_dir<-"/project/getpop/metadata/enviro_data/Envirem_NAmerica/envirem_NON_elev/"

setwd(main_dir)

# read in the coordinates:
coords <- read.csv("/pfs/tc1/project/getpop/metadata/getula_genome_coords.csv")




#### Extract ecological data
###########################################################################

# Pull out altitude and climatic data for coordinates across all species. Do this only once and then write
#    a csv of the data to be used later
# bioclim data is downloaded from: https://www.worldclim.org/data/worldclim21.html  -- Bioclimatic variables at 30s -- specific DL link: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip
#    altitude also from https://www.worldclim.org/data/worldclim21.html -- specific download link: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip


setwd(worldclim_dir)
alt_rast<-rast("wc2.1_30s_elev.tif")  # bring in altitude data
all_alt<-extract(alt_rast, coords[,c("lon", "lat")])[,2] # get the altitiude data for each coordinate for all species
## check if there are any NA values
na_index<-which(is.na(all_alt)) 
na_index # looks good
# get the row indices of NA values
rm(alt_rast) # remove alt_rast object now that data has been extracted from it
names(all_alt) <- coords[,"number"]
setwd(main_dir)
all_alt<-as.matrix(all_alt)
colnames(all_alt)<-"alt"
write.csv(all_alt, file="getula_coords_altitude.csv")

## now do the bioclim layers - current ecological data at 30s resolution
setwd(bioclim_dir)
bioclim_files <- list.files(pattern = "\\.tif$") # list out all files
cur_bioclim_stack<-rast(bioclim_files) # bring the data in as a raster stack
all_cur_bioclim<-extract(cur_bioclim_stack, coords[,c("lon", "lat")])
all_cur_bioclim <- all_cur_bioclim[, !colnames(all_cur_bioclim) == "ID"]
##  Downstream of this, found that there are NA values from samples offshore
colnames(all_cur_bioclim)[colSums(is.na(all_cur_bioclim)) > 0] ## looks like no NA, that's good
rm(cur_bioclim_stack) # remove cur_bioclim_stack from environment now that data has been extracted from it
rownames(all_cur_bioclim) <- coords[,"number"]
setwd(main_dir)
write.csv(all_cur_bioclim, file="getula_coords_present_bioclim.csv")


## Do the same thing for Envirem layers, described here: http://envirem.github.io/
setwd(envirem_dir)
envirem_files <- list.files(pattern = "\\.tif$", recursive=TRUE) # list out all files
cur_envirem_stack <- rast(envirem_files) # bring the data in as a raster stack
all_cur_envirem <- extract(cur_envirem_stack, coords[,c("lon", "lat")])
all_cur_envirem <- all_cur_envirem[, !colnames(all_cur_envirem) %in% "ID"]
##  check for NAs
colnames(all_cur_envirem)[colSums(is.na(all_cur_envirem)) > 0] ## looks good
rm(cur_envirem_stack) # remove cur_envirem_stack from environment now that data has been extracted from it
rownames(all_cur_envirem)<-coords[,"number"]
setwd(main_dir)
write.csv(all_cur_envirem, file="getula_coords_present_envirem.csv")




## Look at correlations 


pdf(file="bioclim_corrs.pdf", width=50, height=50)
pairs.panels(all_cur_bioclim, scale=T) # find correlations
dev.off()


# to drop - drop any correlated above 0.8
bio_to_drop <- c("wc2.1_30s_bio_10", "wc2.1_30s_bio_11", "wc2.1_30s_bio_6", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19", "wc2.1_30s_bio_2", "wc2.1_30s_bio_7")
red_bioclim <- all_cur_bioclim[,!colnames(all_cur_bioclim) %in% bio_to_drop] # reduced set of bioclim variables, still for all species

# combine reduced bioclims with envirems
red_bioc_all_envirem<-cbind(red_bioclim, all_cur_envirem)


## Now see to what degree envirem predictors are highly correlated with these or others
pdf(file="all_envirem_corrs.pdf", width=60, height=60)
pairs.panels(red_bioc_all_envirem, scale=T) # For all bioclims across all points for all species, find correlations
  # do this once and then use just the same set of environmental variables for all - allows more direct comparison
dev.off()

# to drop - drop any correlated above 0.8
envirem_to_drop<-c("current_30arcsec_PETColdestQuarter", "current_30arcsec_PETDriestQuarter", "current_30arcsec_PETWarmestQuarter", "current_30arcsec_PETWettestQuarter", "current_30arcsec_aridityIndexThornthwaite", "current_30arcsec_climaticMoistureIndex", "current_30arcsec_continentality", "current_30arcsec_embergerQ", "current_30arcsec_growingDegDays0", "current_30arcsec_growingDegDays5", "current_30arcsec_maxTempColdest", "current_30arcsec_monthCountByTemp10", "current_30arcsec_thermicityIndex")
bio_envir_red <- red_bioc_all_envirem[,!colnames(red_bioc_all_envirem) %in% envirem_to_drop]


# double check that we don't have any high correlations here
pdf(file="reduced_bio_evirem.pdf", width=60, height=60)
pairs.panels(bio_envir_red, scale=T) # For all bioclims across all points for all species, find correlations
  # do this once and then use just the same set of environmental variables for all - allows more direct comparison
dev.off()
### Looks solid - nothing above 0.8



### Rename a bunch of things in this object
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_1$", "AnnMeanT_B1", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_12", "AnnPrecip_B12", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_15", "PrecipSeas_B15", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_4", "TSeas_B4", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_5", "MaxTWarmMon_B5", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_8", "MTWetQ_B8", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("wc2.1_30s_bio_9", "MTDryQ_B9", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("current_30arcsec_PETseasonality", "PETseasonality", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("current_30arcsec_annualPET", "annualPET", colnames(bio_envir_red))
colnames(bio_envir_red)<-gsub("current_30arcsec_minTempWarmest", "minTempWarmest", colnames(bio_envir_red))

write.csv(bio_envir_red, file="getula_coords_red_enviro_vars.csv")









