## Script to run LFMM on getula genomes and identify any loci that might
##    be under environmental selection

# load up stuff on cluster:
# module load gcc/12.2.0 r/4.4.0 rstudio
# module load udunits/2.2.28 gdal/3.6.1
# make sure base conda isn't activated/doesn't have any weird stuff in it


# load necessary libraries
library(LEA)
library(tidyr)
library(ggplot2)
library(data.table)


# set up some paths
lfmm_dir <- "/pfs/tc1/project/getpop/lfmm/"
# ped_file <- "filtered_ratmapAll_NConly.ped"
env_file <- "/pfs/tc1/project/getpop/metadata/enviro_data/getula_coords_red_enviro_vars.csv"

# ## test on a smaller chromosome at first:
# ped_file <- "filtered_sort_rat_map_all_NC_045557.1_RagTag.ped"

# test on small chromosome with no missing data
# ped_file <- "no_missingNC_045554.ped"

# run on imputed from pca_snmf_getula.R
setwd(lfmm_dir)
ind_names <- readLines("ind_names_lfmm.txt")


# read in LFMM inputed file 
lfmm_file <- "arimap_10k_genome.lfmm_imputed.lfmm"

# get the number of variants for bonferroni correction below
num_variants <- ncol(fread(lfmm_file))

# read in environmental data:
env_data <- read.csv(env_file)

# put enviro data into the order of the genetic data:
env_data <- env_data[match(ind_names, env_data$X),]


# variables to start with: AnnMeanT_B1 (annual mean temp), AnnPrecip_B12 (annual precipitation)

### Precipitation
lfmm_precip_res <- lfmm2(input = lfmm_file, env = env_data$AnnPrecip_B12, K = 3)
# find and plot significant loci
pv <- lfmm2.test(object = lfmm_precip_res,
                 input = lfmm_file,
                 env = env_data$AnnPrecip_B12,
                 full = TRUE)

# Plot
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.05), lty = 2, col = "orange")

# Convert to FDR 
adj_pv <- p.adjust(pv$pvalues, method = "BH")

# Plot out the adjusted p values
plot(-log10(adj_pv), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.10), lty = 2, col = "orange")


## Mean annual temperature
lfmm_temp_res <- lfmm2(input = lfmm_file, env = env_data$AnnMeanT_B1, K = 3)
# find and plot significant loci
pv_temp <- lfmm2.test(object = lfmm_temp_res,
                 input = lfmm_file,
                 env = env_data$AnnMeanT_B1,
                 full = TRUE)

plot(-log10(pv_temp$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.05), lty = 2, col = "orange")


# Convert to FDR 
adj_pv_temp <- p.adjust(pv_temp$pvalues, method = "BH")

# Plot out the adjusted p values
plot(-log10(adj_pv_temp), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.10), lty = 2, col = "orange")






# make a dataframe with p values for each predictor
p_val_df <- data.frame(index = names(adj_pv), p_precip = -log10(adj_pv), p_temp = -log10(adj_pv_temp), row.names = NULL)
p_val_df$index <- rownames(p_val_df)
pval_df_long <- p_val_df %>% pivot_longer(cols=c('p_precip', 'p_temp'),
                    names_to='predictor',
                    values_to='log_p_val')





# ggplot
custom_labels <- c("p_precip" = "Mean annual precipitation",
                   "p_temp" = "Mean annual temperature")


temp_precip_plot <- ggplot(data = pval_df_long, mapping = aes(x = index, y = log_p_val)) +
  geom_point() +
  facet_wrap(~ predictor, ncol = 1, labeller = as_labeller(custom_labels)) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(face="bold"),
        axis.text.y=element_text(size=10)) +
  ylab("-log(10) FDR") +
  xlab("site in genome")

pdf(file = "temp_precip_lfmm.pdf", height = 5, width = 7.5)
print(temp_precip_plot)
dev.off()



png(file = "temp_precip_lfmm.png", height = 5, width = 7.5, units = "in", res = 750)
print(temp_precip_plot)
dev.off()









           