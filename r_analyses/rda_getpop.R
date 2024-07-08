## Script to run RDA on the getula data


# load up stuff on cluster:
# module load gcc/12.2.0 r/4.4.0 rstudio


## load necessary libraries

library(LEA)
library(tidyr)
library(ggplot2)
library(data.table)
library(vegan)


## set up some paths
lfmm_dir <- "/project/getpop/lfmm"
lfmm_file <- "/project/getpop/lfmm/arimap_10k_genome.lfmm_imputed.lfmm"
env_file <- "/pfs/tc1/project/getpop/metadata/enviro_data/getula_coords_red_enviro_vars.csv"
indnames_file <- "ind_names_lfmm.txt"

setwd(lfmm_dir)
ind_names <- readLines(indnames_file)
env_data <- read.csv(env_file)

# put enviro data into the order of the genetic data:
env_data <- env_data[match(ind_names, env_data$X),]

lfmm.data <- read.lfmm(lfmm_file)

# ## remove invariants  - later check if invariants are bc of missing data only ## 
# count_unique <- function(column){
#   length(unique(na.omit(column)))
# } # fxn to count unique characters in a column, excluding NAs
# 
# unique_counts <- apply(lfmm.data, 2, count_unique) # run the fxn on all columns
# 
# columns_tokeep <- unique_counts >= 2 # find which columns have 2 or more unique characters
# lfmm.data_filt <- lfmm.data[, columns_tokeep] #filter based on this
# 
# write.lfmm(lfmm.data_filt, output.file = "testNoinvars.lfmm") # write a new lfmm file
# lfmm_file <- "testNoinvars.lfmm" #reassign it as the new version
# num_variants <- ncol(lfmm.data_filt) # count variants



###############################################
######### Redundancy Analysis (RDA) ###########
###############################################

res_rda <- rda(lfmm.data ~ ., data = env_data[,-1], scale = T)

# adj R^2 
RsquareAdj(res_rda)

summary(res_rda)$concont
screeplot(res_rda)

# check significance for rda  
anova.cca(res_rda, by = "axis")

# variance inflation factors 
vif.cca(cusickii_rda)

plot(cusickii_rda, scaling=3)

load.rda <- summary(cusickii_rda)$species[,1:3]
hist(load.rda[,1], main="Loadings on RDA1", breaks = 20)
hist(load.rda[,2], main="Loadings on RDA2",breaks = 20)
hist(load.rda[,3], main="Loadings on RDA3", breaks = 20) 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) ## 3.
cand2 <- outliers(load.rda[,2], 3) ## 6.
cand3 <- outliers(load.rda[,3], 3) ## 3.

cusickii_rda_cand <- c(names(cand1), names(cand2), names(cand3)) ## j.st the names of the candidates

length(cusickii_rda_cand[duplicated(cusickii_rda_cand)]) ## 7.duplicate detections (detected on multiple RDA axes)
cusickii_rda_cand <- cusickii_rda_cand[!duplicated(cusickii_rda_cand)] ## 1.4 unique candidates 
cusickii_rda_cand

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(lfmm.data_filt) %in% cusickii_rda_cand, 'gray32', '#00000000')
snpcol <- ifelse(colnames(lfmm.data_filt) %in% cusickii_rda_cand, 'red', '#00000000')

## a.es 1 & 2 - zooming in to just the SNPs here...
plot(cusickii_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="cusickii RDA, axes 1 and 2")
points(cusickii_rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(cusickii_rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(cusickii_rda, scaling=3, display="bp", col="#0868ac", cex=1)

## a.es 2 & 3
pdf("RDAaxes23.pdf", height= 8, width = 8)
plot(cusickii_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="cusickii RDA, axes 2 and 3")
points(cusickii_rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(cusickii_rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
text(cusickii_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))
dev.off()

intersetcor(cusickii_rda)[,1:5]


## a.es 4 & 5
plot(cusickii_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(4,5), main="cusickii RDA, axes 4 and 5")
points(cusickii_rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(4,5))
points(cusickii_rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(4,5))
text(cusickii_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(4,5))
