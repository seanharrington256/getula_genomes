## Script to run LFMM on getula genomes and identify any loci that might
##    be under environmental selection

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

# test on a single chrom imputed in pca_snmf_getula.R
setwd(lfmm_dir)
ind_names <- readLines("ind_names_lfmm.txt")


# convert the ped file
# ped2lfmm(ped_file, force = FALSE)
# lfmm_file <- gsub("ped$", "lfmm", ped_file)
lfmm_file <- "rat_map_all_NC_045557.lfmm_imputed.lfmm"

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
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.05/num_variants), lty = 2, col = "orange")


## Mean annual temperature
lfmm_temp_res <- lfmm2(input = lfmm_file, env = env_data$AnnMeanT_B1, K = 3)
# find and plot significant loci
pv_temp <- lfmm2.test(object = lfmm_temp_res,
                 input = lfmm_file,
                 env = env_data$AnnMeanT_B1,
                 full = TRUE)
plot(-log10(pv_temp$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.05/num_variants), lty = 2, col = "orange")


# make a dataframe with p values for each predictor
p_val_df <- data.frame(index = names(pv$pvalues), p_precip = -log10(pv$pvalues), p_temp = -log10(pv_temp$pvalues), row.names = NULL)
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
  geom_hline(yintercept = -log10(0.05/num_variants), linetype = "dashed", color = "red") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(face="bold"),
        axis.text.y=element_text(size=10)) +
  ylab("-log(10) p value") +
  xlab("site on chromosome")

pdf(file = "temp_precip_lfmm.pdf", height = 3, width = 6)
print(temp_precip_plot)
dev.off()








           