### An executable R script to plot pixy output for a bunch of independently analyzed chromosomes
### usage: Rscript plot_pixy_cmd.R <pixy_files_directory> <output_directory>


# Necessary Packages
library(tidyverse)



# print out the command line args used
print(paste0("pixy direectory = ", commandArgs(trailingOnly = TRUE)[[1]]))
print(paste0("output directory = ", commandArgs(trailingOnly = TRUE)[[2]]))



# pixy_dir <- "~/Desktop/pixy_test"
# out_dir <- "~/Desktop/pixy_test/pixy_plots"

# define the input and output directories
pixy_dir <- commandArgs(trailingOnly = TRUE)[[1]]
out_dir <- commandArgs(trailingOnly = TRUE)[[2]]


# Set Working directory
setwd(pixy_dir)



#### pixy conversion function
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
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

#### pixy plot function:
plot_pixy_scaf <- function(pixy_df){
  # split avg_pi statistic into separate pi statistics for each population
  pixy_df[pixy_df$statistic == "avg_pi" & pixy_df$pop1 == pop_names[1],]$statistic <- paste0(pop_names[1], "_pi")
  pixy_df[pixy_df$statistic == "avg_pi" & pixy_df$pop1 == pop_names[2],]$statistic <- paste0(pop_names[2], "_pi")
  
  # make named label vector:
  labels_vec <- c(paste0("pi[", pop_names[1], "]"),
                  paste0("pi[", pop_names[2], "]"),
                  "D[XY]",
                  "F[ST]")
  names(labels_vec) <- c(paste0(pop_names[1], "_pi"), paste0(pop_names[2], "_pi"), "avg_dxy", "avg_wc_fst")
  
  
  ## Plot it?????
  # custom labeller for special characters in pi/dxy/fst
  pixy_labeller <- as_labeller(labels_vec,
                               default = label_parsed)
  
  if(length(unique(pixy_df$chromosome)) > 1){
    print("MORE THAN ONE CHOMROSOME!!!!")
  }else{
    scaf_name <- pixy_df$chromosome[[1]]
    
    # plotting summary statistics along a single chromosome
    pixy_df %>%
      # filter(chromosome == "NC_045554.1_RagTag") %>%
      filter(statistic %in% c(paste0(pop_names[1], "_pi"), paste0(pop_names[2], "_pi"), "avg_dxy", "avg_wc_fst")) %>%
      mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
      ggplot(aes(x = chr_position, y = value, color = statistic))+
      geom_line(size = 0.25)+
      facet_grid(statistic ~ .,
                 scales = "free_y", switch = "x", space = "free_x",
                 labeller = labeller(statistic = pixy_labeller,
                                     value = label_value))+
      xlab(paste0("Position on ", scaf_name, " (Mb)"))+
      ylab("Statistic Value")+
      theme_bw()+
      theme(panel.spacing = unit(0.1, "cm"),
            strip.background = element_blank(),
            strip.placement = "outside",
            legend.position = "none")+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      scale_color_brewer(palette = "Set1")
  }
}




# read in pixy files and generate a dataframe
pixy_files <- list.files(pattern = "_pi\\.txt$|_dxy\\.txt$|_fst\\.txt$", full.names = TRUE)


# ADDED - Extract unique chromosome names from the filenames
chromosome_names <- unique(gsub("pi\\.txt$|dxy\\.txt$|fst\\.txt$", "", basename(pixy_files)))
# Extract filenames "basename" 
# "gsub" used for pattern matching and replacement 
# "Unique only gets unique elements from the resulting list



## loop over chromosomes/scaffolds and plot to pdf named for current directory
current_dir <- basename(pixy_dir)

pdf(file = paste0(out_dir, "/", current_dir, ".pdf"), width = 8, height = 6.5)
for(scaf in chromosome_names){ # for each scaffold in the unique chromosomes/scaffolds
  
  scaf_files <- grep(scaf, pixy_files, value = TRUE) # get the files for the current scaffold
  pixy_df <- pixy_to_long(scaf_files) # read in the current scaffold
  pop_names <- unique(pixy_df$pop1) # # get names of the populations
  print(plot_pixy_scaf(pixy_df))
  
}
dev.off()
