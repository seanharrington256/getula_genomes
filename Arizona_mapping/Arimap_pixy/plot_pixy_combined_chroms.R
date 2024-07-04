# On beartooth, need to run:
# module load gcc/12.2.0 r/4.4.0 rstudio


# Necessary Packages
library(tidyverse)

# Set Working directory

## Set up file paths for 100 kb windows
# pixy_dir <- "/pfs/tc1/project/getpop/pixy_out/pixy_out100_combined/"
# out_dir <- "/pfs/tc1/project/getpop/pixy_out/plots_pixy100_combined/"
# outfile_name <- "pixy100_combined"

## Set up file paths for 10 kb windows
pixy_dir <- "/project/getpop/pixy_out_AriMap/pixy_out10_combined"
out_dir <- "/project/getpop/pixy_out_AriMap/pixy_out10_combined"
outfile_name <- "pixy10_combined"


setwd(pixy_dir)



# # Define the function to extract unique chromosome names
# extract_chromosome_names <- function(pixy_files) {
#   unique(gsub("pi\\.txt$|dxy\\.txt$|fst\\.txt$", "", basename(pixy_files)))
# }

# Define the function to convert pixy files to long format
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  
  for (i in 1:length(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if (stat_file_type == "pi") {
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
    } else {
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

# # Function to plot summary statistics across all chromosomes
# plot_genome_wide_pixy <- function(pixy_df) {
#   # Plotting summary statistics across all chromosomes
#   pixy_df %>%
#     mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
#                                          chromosome == "X" ~ "even",
#                                          TRUE ~ "odd" )) %>%
#     mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
#     filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
#     ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
#     geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
#     facet_grid(statistic ~ chromosome,
#                scales = "free_y", switch = "x", space = "free_x",
#                labeller = labeller(statistic = pixy_labeller,
#                                    value = label_value))+
#     xlab("Chromsome")+
#     ylab("Statistic Value")+
#     scale_color_manual(values = c("grey50", "black"))+
#     theme_classic()+
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.spacing = unit(0.1, "cm"),
#           strip.background = element_blank(),
#           strip.placement = "outside",
#           legend.position ="none")+
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
# }


# Custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi[within]",
                               avg_dxy = "pi[between]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)



# Read in pixy files and generate a dataframe
pixy_files <- list.files(pattern = "_pi\\.txt$|_dxy\\.txt$|_fst\\.txt$", full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

# Make a new column for chromosome that is just a number
pixy_df$chrom_name <- pixy_df$chromosome
pixy_df$chromosome <- as.numeric(as.factor(pixy_df$chromosome))

# Let's get this down to just 18 chromosomes:
pixy_df <- filter(pixy_df, chromosome < 19)


# get out Fst values and find upper quantile
all_fst <- filter(pixy_df, statistic == "avg_wc_fst")$value
upper_99th <- quantile(all_fst, probs = 0.99, na.rm = TRUE)

# make a column that is chromosome and position together
pixy_df <- mutate(pixy_df, chr_pos = paste(chromosome, window_pos_1, sep = "_"))

# get chr_pos of top 1% Fst outliers
top_01_fs_twinds <- filter(pixy_df, statistic == "avg_wc_fst" & value >= upper_99th)$chr_pos

# Find upper 1% dxy outliers
all_dxy <- filter(pixy_df, statistic == "avg_dxy")$value
upper_99dxy <- quantile(all_dxy, probs = 0.99, na.rm = TRUE)


# get chr_pos of top 1% Fst AND dxy outliers
top_01_fst_dxy_winds <- filter(pixy_df, statistic == "avg_dxy" & value >= upper_99dxy & chr_pos %in% top_01_fs_twinds)$chr_pos


# ## make the outliers red
# plot_genome_wide_pixy <- function(pixy_df, outlier_winds) {
#   # Plotting summary statistics across all chromosomes
#   pixy_df %>%
#     mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
#                                          chromosome == "X" ~ "even",
#                                          TRUE ~ "odd" )) %>%
#     
#     mutate(chrom_color_group = ifelse(chr_pos %in% outlier_winds, "outlier", chrom_color_group)) %>%
#     mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
#     filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
#     ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
#     geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
#     facet_grid(statistic ~ chromosome,
#                scales = "free_y", switch = "x", space = "free_x",
#                labeller = labeller(statistic = pixy_labeller,
#                                    value = label_value))+
#     xlab("Chromsome")+
#     ylab("Statistic Value")+
#     scale_color_manual(values = c("grey50", "black", "red"))+ # Adding "red" color for outliers
#     theme_classic()+
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.spacing = unit(0.1, "cm"),
#           strip.background = element_blank(),
#           strip.placement = "outside",
#           legend.position ="none")+
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
# }





# ## make the outliers red and large
# plot_genome_wide_pixy <- function(pixy_df, outlier_winds){
#   pixy_df %>%
#     mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
#                                          chromosome == "X" ~ "even",
#                                          TRUE ~ "odd" )) %>%
#     # mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
#     mutate(chromosome = factor(chromosome, levels = 1:18)) %>% ## Editing this for snake chromosome number - excluding small little chrs
#     mutate(chrom_color_group = ifelse(chr_pos %in% outlier_winds, "outlier", chrom_color_group)) %>%
#     filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
#     ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
#     geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
#     geom_point(data = . %>% filter(chrom_color_group == "outlier"), # Filtered for outliers
#                size = 1, # Set size for outliers
#                color = "red", alpha = 0.75, stroke = 0) + # Set color for outliers
#     facet_grid(statistic ~ chromosome,
#                scales = "free", switch = "x", space = "free_x",
#                labeller = labeller(statistic = pixy_labeller,
#                                    value = label_value))+
#     xlab("Chromosome")+
#     ylab("Statistic Value")+
#     scale_color_manual(values = c("grey50", "black", "red"))+ # Adding "red" color for outliers
#     theme_classic()+
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.spacing = unit(0.1, "cm"),
#           strip.background = element_blank(),
#           strip.placement = "outside",
#           legend.position ="none",
#           strip.text.x = element_text(size = 6))+
#     scale_x_continuous(expand =  expansion(add = 10000000)) +
#     scale_y_continuous(expand = c(0, 0.005), limits = c(0, NA))
# }


## make the outliers red and large - blue if outlier for Fst & dxy
plot_genome_wide_pixy <- function(pixy_df, outlier_winds, fst_dxy_outwinds){ # takes a pixy dataframe, a vector of windows that are fst outliers, and a vector of windows that are both Fst AND dxy outliers
  pixy_df %>%
    mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                         chromosome == "X" ~ "even",
                                         TRUE ~ "odd" )) %>%
    # mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
    mutate(chromosome = factor(chromosome, levels = 1:18)) %>% ## Editing this for snake chromosome number - excluding small little scaffolds
    mutate(chrom_color_group = ifelse(chr_pos %in% outlier_winds, "outlier", chrom_color_group)) %>% # change color groups for Fst outliers
    mutate(chrom_color_group = ifelse(chr_pos %in% fst_dxy_outwinds, "fstdxyoutlier", chrom_color_group)) %>%  # change color groups for Fst AND dxy outliers
    filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
    ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
    geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
    geom_point(data = . %>% filter(chrom_color_group == "outlier"), # Filtered for outliers
               size = 1, # Set size for outliers
               color = "forestgreen", alpha = 0.75, stroke = 0) + # Set color for outliers
    geom_point(data = . %>% filter(chrom_color_group == "fstdxyoutlier"), # Filtered for Fst AND dxy outliers
               size = 3, # Set size for outliers
               color = "blue", alpha = 0.8, stroke = 0) + # Set color for fst & dxy outliers
    facet_grid(statistic ~ chromosome,
               scales = "free", switch = "x", space = "free_x",
               labeller = labeller(statistic = pixy_labeller,
                                   value = label_value))+
    xlab("Chromosome")+
    ylab("Statistic Value")+
    scale_color_manual(values = c("grey50", "forestgreen", "black", "blue"))+ # Adding "red" and blue color for outliers
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position ="none",
          strip.text.x = element_text(size = 6))+
    scale_x_continuous(expand =  expansion(add = 10000000)) +
    scale_y_continuous(expand = c(0, 0.005), limits = c(0, NA))
}


# write to file
outfile <- paste0(out_dir, outfile_name, ".pdf")
pdf(file = outfile, width = 8, height = 6.4)
plot_genome_wide_pixy(pixy_df, top_01_fs_twinds, top_01_fst_dxy_winds)
dev.off()

#plot png
outfile_png <- paste0(out_dir, outfile_name, ".png")
png(file = outfile_png, width = 8, height = 6.5, units = "in", res = 750)
plot_genome_wide_pixy(pixy_df, top_01_fs_twinds, top_01_fst_dxy_winds)
dev.off()

