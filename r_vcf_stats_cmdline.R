#!/bin/bash R

#View in R
##load tidyverse
library(tidyverse)

# script to plot out summaries of vcf file
# runs from the command line, and takes 2 arguments
# first arg is the direcgory containing the summary files
# second arg is the basename of summary files


vcf_stats_dir <- commandArgs(trailingOnly = TRUE)[[1]] 
file_basename <- commandArgs(trailingOnly = TRUE)[[2]] 


setwd(vcf_stats_dir)




#view stats
var_qual <- read_delim(paste0(file_basename, ".lqual"), delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)
#plot data
variant_qual <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

rm(var_qual)

#variant mean depth
var_depth <- read_delim(paste0(file_basename, ".ldepth.mean"), delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

#plot data
variant_depth <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + xlim(0, 100)

#Stats
var_dep <- summary(var_depth$mean_depth)

rm(var_depth)


#variant missingness
var_miss <- read_delim(paste0(file_basename, ".lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
#plot
variant_miss<- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#stats
miss <- summary(var_miss$fmiss)

rm(var_miss)

#allele frequency for minor
var_freq <- read_delim(paste0(file_basename, ".frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
                       
##minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

#plot
variant_freq <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#view data
frequency <- summary(var_freq$maf)

rm(var_freq)

#individual based stats


##mean depth/individual
ind_depth <- read_delim(paste0(file_basename, ".idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
#plot
indiv_depth <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)


#proportion of missing data/individual
ind_miss  <- read_delim(paste0(file_basename, ".imiss"), delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
#plot
indiv_miss<- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)


#heterozygosity and inbreeding coef/individual
ind_het <- read_delim(paste0(file_basename, ".het"), delim = "\t",
           col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
           
#plot
indiv_het <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)

paste0("vcf_plots_", file_basename, ".pdf")

pdf(file=paste0(paste0(file_basename, "_vcf_plots.pdf")), width=8, height=5)
plot(variant_qual, col="lightblue", cex=1.2, pch = 19)
plot(variant_depth, col="lightblue", cex=1.2, pch = 19)
plot(variant_freq, col="lightblue", cex=1.2, pch = 19)
plot(indiv_depth, col="lightblue", cex=1.2, pch = 19)
plot(indiv_miss, col="lightblue", cex=1.2, pch = 19)
plot(indiv_het, col="lightblue", cex=1.2, pch = 19)
dev.off()



combined_stats <- rbind(var_dep, miss, frequency)
write.csv(combined_stats, file=paste0("vcf_stats_", file_basename, ".csv"))

