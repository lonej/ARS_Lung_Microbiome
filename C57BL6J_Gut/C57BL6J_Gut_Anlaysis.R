#!/usr/bin/env Rscript

# --------------------------------------------------
# Title : Main_Script - C57BL6J GUT Analysis
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Script
# --------------------------------------------------
# Objective : Run all functions
# --------------------------------------------------
#
# Main script
#
# --------------------------------------------------




# ##################################################
# Load Packages
# ##################################################

library(tidyverse)
library(pheatmap)
library(vegan)
library(ggcorrplot)
library(ggplot2)
library(lefser)
library(SummarizedExperiment)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(fossil)
library(skimr)
library(reshape2)
library(ggpubr)
library(svglite)
library(factoextra)




# ##################################################
# Directories structures check
# ##################################################

# Check the existence of all needed directories
for (dir in c("Checkpoints_Cleaning",
              "Output_Adonis",
              "Rarefied_Dataset")){
  
  # Delete directories (and files) if they exist
  # unlink(paste0("./", dir, "/", sep=""), recursive = TRUE)
  
  # Create the directories if it does not exist
  if (!(file.exists(paste0("./", dir, "/", sep="")))){
    dir.create(paste0("./", dir, "/", sep=""))
    print(paste(dir, "Created", sep = " ->  "))
  } else if (file.exists(paste0("./", dir, "/", sep=""))){
    print(paste(dir, "Ready", sep = " ->  "))
  }
  
}







# ##################################################
# Rarefaction and data cleaning
# ##################################################

# Preparation
rm(list=ls())
source("./Functions/Clean_N_Rarefied.R")
set.seed(4242)

# Parameters
min_reads <- 20000
rare_curve <- FALSE # TRUE or FALSE

# Cleaning datasets and Rarefy
Clean_N_Rarefied(min_reads, rare_curve)







# ##################################################
# PERMANOVA of rarefied Dataset
# ##################################################

# Preparation
rm(list=ls())
source("./Functions/Adonis_Loop.R")

# Set seed
set.seed(4242)

# Options for running
distance_ecology <- "bray"
graph_or_not <- "Yes"
taxonomy_rank <- "OTU"
medium <- "IMM"

# Load Rarefied Dataset
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
                                 sep=""), check.names=FALSE)

# Removing Informations (Day 12 and Microbiome colums)
exp_OTU_Table <- exp_OTU_Table %>% filter(day_postInfection !=12)


for (genotype in c("C57BL6J")){
  
  # Initialize Data frame for loop 
  Data <- data.frame(Mouse_strain=NULL, Culture=NULL, Days=NULL, 
                     Compare=NULL, Pval = NULL, R2 = NULL, Betadisp = NULL)
  
  # List of comparison for adonis
  if (genotype == "C57BL6J"){
    parasite_list <- c("NA", "NK", "AK", "NAK")
  } else if (genotype == "DBA"){
    parasite_list <- c("NA", "NS", "AS", "NAS")
  }
  
  # Run PERMANOVA on all days and comparison list
  for (dpi in c(5)){
    for (pair_parasite in parasite_list){
      res <- Adonis_Loop(data=exp_OTU_Table, culture=medium, Days=dpi,
                         mouse=genotype, parasite=pair_parasite,
                         taxo_rank=taxonomy_rank, distance=distance_ecology, 
                         graph=graph_or_not)
      Data <- rbind.data.frame(Data, res)
    }
  }
  
  # Round R2 and Betadisp and reorder according to what we generated previously
  Data <- Data %>% mutate(across(6:7, round, 3)) %>%
    arrange(Culture, factor(Days, levels = c(5)))
  
  # Save the Datset
  write.csv2(Data, paste("./Output_Adonis/", 
                         paste(taxonomy_rank, genotype, medium, distance_ecology, "Pval", sep="_"),
                         ".csv", sep=""), row.names = TRUE)
}






