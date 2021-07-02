#!/usr/bin/env Rscript

# --------------------------------------------------
# Title : Main_Script - DBA Lung Analysis
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Script
# --------------------------------------------------
# Objective : Run all functions
# --------------------------------------------------
#
# Main Script
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
              "Checkpoints_Sentinels_and_Blanck",
              "Output_Adonis",
              "Output_Diversity",
              "Output_Heatmaps",
              "Output_LEfSe",
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
min_reads <- 5000
rare_curve <- FALSE # Compute Rarefaction curve TRUE or FALSE

# Cleaning datasets and Rarefy
Clean_N_Rarefied(min_reads, rare_curve)





# ##################################################
# Find outlier Sentinels
# ##################################################
# 
# # Preparation
rm(list=ls())
source("./Functions/Sentinel_Outliers.R")

# Parameters
taxonomy_rank <- "OTU"

# Load Rarefied Dataset
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
                                 sep=""), check.names=FALSE)
exp_OTU_Table <-  exp_OTU_Table %>% dplyr::select(-Microbiome)

# List of outlier cages to remove
Sentinel_Outliers(data = exp_OTU_Table,
                  taxonomy_rank = taxonomy_rank)






# ##################################################
# Position sample to Blanck
# ##################################################

# Preparation
rm(list=ls())
source("./Functions/Blanck_PCoA.R")

# Parameters
taxonomy_rank <- "OTU"
mouse2rm <- c("D7_2", "D7_7" , "D7_10")

# Load Dataset with Blanck
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_WithBlanck_",taxonomy_rank,"_Table.csv",
                                 sep=""), check.names=FALSE)
data <-  exp_OTU_Table %>% dplyr::select(-Microbiome)

# PCoA to see if Blanck is different
Blanck_PCoA(data, mouse2rm)







# ##################################################
# Convert ANKA to ANKA_lo and ANKA_hi
# ##################################################

# Clean 
rm(list=ls())

for (taxonomy_rank in c("domain", "phylum", "class",
                        "order", "family", "genus", "OTU")) {
  
  exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
                                   sep=""), check.names=FALSE)
  parasitaemia <- read.csv2("./Dataset/Parasitaemia.csv")
  
  # Find common value and remove the filtered data from parasitaemia
  parasitaemia <- parasitaemia %>% 
    filter(sample.id %in% intersect(parasitaemia$sample.id, exp_OTU_Table$ID) )
  
  # Rename the ANKA with the lo or hi
  exp_OTU_Table$parasite_strain[exp_OTU_Table$ID %in% parasitaemia$sample.id] <- parasitaemia$Group 
  
  # Save it!
  write.csv2(exp_OTU_Table, paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,
                                  "_Table.csv", sep=""),row.names=F)
}





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
genotype <-  "DBA"
rm_list <- c("D7_2", "D7_7" , "D7_10")

# Load Rarefied Dataset
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
                                 sep=""), check.names=FALSE)

# Removing outliers
exp_OTU_Table <- exp_OTU_Table %>% 
  dplyr::filter(!(mouse_cage %in% rm_list)) %>% 
  dplyr::select(-Microbiome)

# Initialize Data frame for loop 
Data <- data.frame(Mouse_strain=NULL, Culture=NULL, Days=NULL, 
                   Compare=NULL, Pval = NULL, R2 = NULL, Betadisp = NULL)

# List of comparison for adonis
parasite_list <- c("NAl", "NAh", "NS", "AlS", "AhS", "AlAh", "NAlAhS")


# Run PERMANOVA on all days and comparison list
for (dpi in c(5, 7)){
  for (pair_parasite in parasite_list){
    
    table1 <- exp_OTU_Table
    res <- Adonis_Loop(data=table1, culture=medium, Days=dpi,
                       mouse=genotype, parasite=pair_parasite,
                       taxo_rank=taxonomy_rank, distance=distance_ecology, 
                       graph=graph_or_not)
    Data <- rbind.data.frame(Data, res$stats_adonis)
  }
}

# Round R2 and Betadisp and reorder according to what we generated previously
Data <- Data %>% mutate(across(6:7, round, 3)) %>%
  arrange(Culture, factor(Days, levels = c(357, 5, 7)))

# Save the Datset
write.csv2(Data, paste("./Output_Adonis/", 
                       paste(taxonomy_rank, genotype, medium, distance_ecology, "Pval", sep="_"),
                       ".csv", sep=""), row.names = TRUE)





# ##################################################
# LEfSe of rarefied Dataset
# ##################################################

# Clear memory
rm(list=ls())
source("./Functions/lefseByPair.R")
rm_list <- c("D7_2", "D7_7" , "D7_10")

# Parameters
taxonomy_rank <- "family"

# Load Dataset
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,
                                 "_Table.csv", sep=""), check.names=FALSE)
exp_OTU_Table <- exp_OTU_Table %>% dplyr::filter(!(mouse_cage %in% rm_list))

# Day 5
exp_OTU_Table5 <- exp_OTU_Table %>% dplyr::select(-c(mouse_cage, Microbiome) ) 
res5 <- lefseByPair(dataset=exp_OTU_Table5, taxo_rank=taxonomy_rank, 
                   mouse="DBA", Days=5)
# Day 7
exp_OTU_Table7 <- exp_OTU_Table %>% dplyr::select(-c(mouse_cage, Microbiome) ) 
res7 <- lefseByPair(dataset=exp_OTU_Table7, taxo_rank=taxonomy_rank, 
                   mouse="DBA", Days=7)









# ##################################################
# Heatmap trimmed of low abundant taxa
# ##################################################


# Reset Memory
rm(list=ls())
source("./Functions/my_Heatmap2.R")
rm_list <- c("D7_2", "D7_7" , "D7_10")

# Parameters
taxonomy_rank <- "family"
genotype <- "DBA"

# Load Datasets
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,
                                 "_Table.csv", sep=""), check.names=FALSE)
exp_OTU_Table <- exp_OTU_Table %>% dplyr::filter(!(mouse_cage %in% rm_list))

# Information selection
exp_OTU_Table <- exp_OTU_Table %>% dplyr::select(-c(Microbiome, mouse_cage))
exp_OTU_Table <- exp_OTU_Table %>% filter(parasite_strain !="Sentinel")

# Do Heatmaps
nrow(exp_OTU_Table %>% dplyr::filter(day_postInfection==5))
my_Heatmap2(data = exp_OTU_Table,
            days=5, max_zero_taxa=22, 
            taxo_rank=taxonomy_rank, mouse=genotype)

nrow(exp_OTU_Table %>% dplyr::filter(day_postInfection==7))
my_Heatmap2(data = exp_OTU_Table,
            days=7, max_zero_taxa=18, 
            taxo_rank=taxonomy_rank, mouse=genotype)


dev.off()





# ##################################################
# Alpha diversity 
# ##################################################

# Clear memory
rm(list=ls())
source("./Functions/my_Diversity.R")
rm_list <- c("D7_2", "D7_7" , "D7_10")

# Parameters
taxonomy_rank <- "OTU"

# Load Dataset
exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
                                 sep=""), check.names=FALSE)
exp_OTU_Table <- exp_OTU_Table %>% dplyr::filter(!(mouse_cage %in% rm_list))
exp_OTU_Table <- exp_OTU_Table %>% 
  dplyr::select(-c(Microbiome)) %>% 
  filter(parasite_strain != "Sentinel")

# Diversity Graphs
for (dpi in c(5, 7)) {
  my_Diversity(table=exp_OTU_Table, taxo_rank=taxonomy_rank, 
               Days=dpi, mouse="DBA", adj_meth = "BH")
}



