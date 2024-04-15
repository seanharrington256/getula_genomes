# Necessary Packages
library(tidyverse)

# Set Working directory
pixy_dir <- "/pfs/tc1/project/getpop/pixy_out/pixy_out100_combined/"
out_dir <- "/pfs/tc1/project/getpop/pixy_out/plots_pixy100_combined/"
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
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)



# Read in pixy files and generate a dataframe
pixy_files <- list.files(pattern = "_pi\\.txt$|_dxy\\.txt$|_fst\\.txt$", full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

# Make a new column for chromosome that is just a number
pixy_df$chrom_name <- pixy_df$chromosome
pixy_df$chromosome <- as.numeric(as.factor(pixy_df$chromosome))

# get out Fst values and find upper quantile
all_fst <- filter(pixy_df, statistic == "avg_wc_fst")$value
upper_99th <- quantile(all_fst, probs = 0.99)

# make a column that is chromosome and position together
pixy_df <- mutate(pixy_df, chr_pos = paste(chromosome, window_pos_1, sep = "_"))

# get chr_pos of top 1% Fst outliers
top_01_fs_twinds <- filter(pixy_df, statistic == "avg_wc_fst" & value >= upper_99th)$chr_pos


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


## make the outliers red and large
plot_genome_wide_pixy <- function(pixy_df, outlier_winds){
  pixy_df %>%
    mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                         chromosome == "X" ~ "even",
                                         TRUE ~ "odd" )) %>%
    mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"))) %>%
    mutate(chrom_color_group = ifelse(chr_pos %in% outlier_winds, "outlier", chrom_color_group)) %>%
    filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
    ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
    geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
    geom_point(data = . %>% filter(chrom_color_group == "outlier"), # Filtered for outliers
               size = 1, # Set size for outliers
               color = "red", alpha = 0.75, stroke = 0) + # Set color for outliers
    facet_grid(statistic ~ chromosome,
               scales = "free", switch = "x", space = "free_x",
               labeller = labeller(statistic = pixy_labeller,
                                   value = label_value))+
    xlab("Chromosome")+
    ylab("Statistic Value")+
    scale_color_manual(values = c("grey50", "black", "red"))+ # Adding "red" color for outliers
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
outfile <- paste0(out_dir, "pixy100_combined.pdf")
pdf(file = outfile, width = 8, height = 6.4)
plot_genome_wide_pixy(pixy_df, top_01_fs_twinds)
dev.off()

#plot png
outfile_png <- paste0(out_dir, "pixy100_combined.png")
png(file = outfile_png, width = 8, height = 6.5, units = "in", res = 750)
plot_genome_wide_pixy(pixy_df, top_01_fs_twinds)
dev.off()

