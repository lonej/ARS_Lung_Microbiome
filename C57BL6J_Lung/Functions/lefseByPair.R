

# --------------------------------------------------
# Title : lefseByPair
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Linear discriminant analysis Effect Size
# --------------------------------------------------
# Input(s) : 
# dataset : OTU table with metadata
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
# mouse : "C57BL6J" or "DBA
# Days : 3, 5 or 357 (for day 5 and 7 combined)
# --------------------------------------------------
# Output(s) : 
# LEfSe plots : Raster and vector graphics (LEfSe results LDA)
# --------------------------------------------------


lefseByPair <- function(dataset, taxo_rank, mouse, Days){
  
  
  # Debugging tools
  # rm(list=ls())
  # taxo_rank <- "family"
  # mouse <- "C57BL6J"
  # Days <- 5
  # dataset <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # dataset <- dataset %>% dplyr::select(-Microbiome)

  
  
  
  # --------------------------------------------------
  # Data frame for pairwise parasite comparison
  # --------------------------------------------------
  
  # Parameters
  LDA_tresh <- 2
  
  # Select the genetic background
  data <- dataset %>% filter(mouse_strain==mouse)
  data <- data %>% column_to_rownames("ID")
  
  # Sub-group construction for BL6
  if (mouse == "C57BL6J") {
    
    NK_group <- data %>% filter(parasite_strain=="NI" |
                                  parasite_strain=="K173",
                                day_postInfection == Days)
    AK_group <- data %>% filter(parasite_strain=="ANKA" |
                                  parasite_strain=="K173",
                                day_postInfection == Days)
    NA_group <- data %>% filter(parasite_strain=="NI" |
                                  parasite_strain=="ANKA",
                                day_postInfection == Days)

  } 
  
  # Re-ordering Factors to compare control (left) not controls (right)
  if (mouse == "C57BL6J") {
    NK_group$parasite_strain <- factor(NK_group$parasite_strain, 
                                       levels = c("NI", "K173"))
    AK_group$parasite_strain <- factor(AK_group$parasite_strain, 
                                       levels = c("ANKA", "K173"))
    NA_group$parasite_strain <- factor(NA_group$parasite_strain, 
                                       levels = c("NI", "ANKA"))
  }
  
  
  
  # --------------------------------------------------
  # Creation of SummarizedExperiment obbject
  # --------------------------------------------------
  
  # Get the Metadata
  meta_NK_group <- NK_group %>% select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  meta_AK_group <- AK_group %>% select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  meta_NA_group <- NA_group %>% select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  
  # Get the count matri
  NK_group <- NK_group %>% select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  AK_group <- AK_group %>% select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  NA_group <- NA_group %>% select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  
  
  # SummarizedExperiment Dataset 
  se_NK_group <- SummarizedExperiment(t(NK_group), colData = meta_NK_group)
  se_AK_group <- SummarizedExperiment(t(AK_group), colData = meta_AK_group)
  se_NA_group <- SummarizedExperiment(t(NA_group), colData = meta_NA_group)
  
  
  
  
  # --------------------------------------------------
  # LEfSe for each mouse strain
  # --------------------------------------------------
  
  # Function to draw a plot usign ggplot2
  source("./Functions/my_lefserPlot.R")
  
  
  if (mouse == "C57BL6J") {
    
    # LEfSe
    res_se_NK_group <- lefser(se_NK_group, groupCol = "parasite_strain", 
                              blockCol = "day_postInfection", lda.threshold = LDA_tresh)
    res_se_AK_group <- lefser(se_AK_group, groupCol = "parasite_strain", 
                              blockCol = "day_postInfection", lda.threshold = LDA_tresh)
    group_colors <- list(NI="white", ANKA="snow3", K173="firebrick1")
    
    
    
    # LEfSE plot for NI vs K173
    plot_NK <- my_lefserPlot(res_se_NK_group, colors = group_colors,
                             meta_levels=levels(meta_NK_group$parasite_strain)) +
      ggtitle("NI vs K173") + 
      theme(plot.title = element_text(size = 20, face = "bold"))
    
    # Save the ggplots !
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"),
                 "_NI_vs_K173.png"), width = 10, height = 6, plot = plot_NK)
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"),
                 "_NI_vs_K173.svg"), width = 10, height = 6, plot = plot_NK)
    
    
    
    # LEfSE plot for ANKA vs K173
    plot_AK <- my_lefserPlot(res_se_AK_group, colors = group_colors, 
                             meta_levels=levels(meta_AK_group$parasite_strain)) +
      ggtitle("ANKA vs K173") + 
      theme(plot.title = element_text(size = 20, face = "bold"))
    
    # Save the ggplots !
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_ANKA_vs_K173.png"), width = 10, height = 6, plot = plot_AK)
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_ANKA_vs_K173.svg"), width = 10, height = 6, plot = plot_AK)
    
    
    
    
  } # Enf if B6
  
} # End of function




