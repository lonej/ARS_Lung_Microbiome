

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
  # mouse <- "DBA"
  # Days <- 5
  # rm_list <- c("D7_2", "D7_7" , "D7_10")
  # dataset <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # dataset <- dataset %>% dplyr::filter(!(mouse_cage %in% rm_list))
  # dataset <- dataset %>% dplyr::select(-c(mouse_cage, Microbiome) )

  
  
  # --------------------------------------------------
  # Data frame for pairwise parasite comparison
  # --------------------------------------------------
  
  # Select the genetic background
  data <- dataset %>% filter(mouse_strain==mouse)
  data <- data %>% column_to_rownames("ID")
  
  # Sub-group construction for DBA
 if (mouse == "DBA"){
    
    AhK_group <- data %>% filter(parasite_strain=="ANKA_hi" |
                                   parasite_strain=="SMAC",
                                 day_postInfection == Days)
    NAh_group <- data %>% filter(parasite_strain=="NI" |
                                  parasite_strain=="ANKA_hi",
                                day_postInfection == Days)
    AlAh_group <- data %>% filter(parasite_strain=="ANKA_lo" |
                                   parasite_strain=="ANKA_hi",
                                 day_postInfection == Days)
  }
  
  # Re-ordering Factors to compare control (left) not controls (right)
 if (mouse == "DBA"){
    AhK_group$parasite_strain <- factor(AhK_group$parasite_strain, 
                                       levels = c("SMAC", "ANKA_hi"))
    NAh_group$parasite_strain <- factor(NAh_group$parasite_strain, 
                                       levels = c("NI", "ANKA_hi"))
    AlAh_group$parasite_strain <- factor(AlAh_group$parasite_strain, 
                                        levels = c("ANKA_lo", "ANKA_hi"))
  }
  
  
  
  # --------------------------------------------------
  # Creation of SummarizedExperiment obbject
  # --------------------------------------------------
  
  # Get the Metadata
  meta_AhK_group <- AhK_group %>% dplyr::select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  meta_NAh_group <- NAh_group %>% dplyr::select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  meta_AlAh_group <- AlAh_group %>% dplyr::select(parasite_strain, mouse_strain, day_postInfection, culture_type)
  
  
  # Get the count matrix
  AhK_group <- AhK_group %>% dplyr::select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  NAh_group <- NAh_group %>% dplyr::select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  AlAh_group <- AlAh_group %>% dplyr::select(-c(parasite_strain, mouse_strain, day_postInfection, culture_type))
  
  
  # SummarizedExperiment Dataset 
  se_AhK_group <- SummarizedExperiment(t(AhK_group), colData = meta_AhK_group)
  se_NAh_group <- SummarizedExperiment(t(NAh_group), colData = meta_NAh_group)
  se_AlAh_group <- SummarizedExperiment(t(AlAh_group), colData = meta_AlAh_group)
  
  
  
  
  # --------------------------------------------------
  # LEfSe for each mouse strain
  # --------------------------------------------------
  
  # Function to draw a plot usign ggplot2
  source("./Functions/my_lefserPlot.R")
  
  if (mouse == "DBA"){
    
    LDA_thresh <- 2
    
    # LEfSe
    res_se_NAh_group <- lefser(se_NAh_group, groupCol = "parasite_strain", 
                               # blockCol = "day_postInfection", 
                               lda.threshold = LDA_thresh)
    res_se_AhK_group <- lefser(se_AhK_group, groupCol = "parasite_strain", 
                               # blockCol = "day_postInfection", 
                               lda.threshold = LDA_thresh)
    res_se_AlAh_group <- lefser(se_AlAh_group, groupCol = "parasite_strain", 
                                # blockCol = "day_postInfection", 
                                lda.threshold = LDA_thresh)
    
    
    
    # Set colors
    group_colors <- list(NI="white", ANKA_lo = "orange", 
                         ANKA_hi ="firebrick3", SMAC="grey70")
    
    
    # LEfSe for NI vs ANKA_hi
    # ------------------------------
    
    # LEfSE plot
    plot_NAh <- my_lefserPlot(res_se_NAh_group, colors = group_colors, 
                             meta_levels=levels(meta_NAh_group$parasite_strain)) +
      ggtitle("NI vs ANKA_hi") + 
      theme(plot.title = element_text(size = 20, face = "bold"))
    plot_NAh
    
    # Save plot 
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_NI_vs_ANKA_hi.png"), width = 10, height = 6, plot = plot_NAh)
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_NI_vs_ANKA_hi.svg"), width = 10, height = 6, plot = plot_NAh)
    
    
    
    # LEfSe for SMAC vs ANKA_hi
    # ------------------------------
    
    # LEfSE plot
    plot_AhK <- my_lefserPlot(res_se_AhK_group, colors = group_colors, 
                              meta_levels=levels(meta_AhK_group$parasite_strain)) +
      ggtitle("SMAC vs ANKA_hi") + 
      theme(plot.title = element_text(size = 20, face = "bold"))
    plot_AhK
    
    # Save plot 
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_SMAC_vs_ANKA_hi.png"), width = 10, height = 6, plot = plot_AhK)
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_SMAC_vs_ANKA_hi.svg"), width = 10, height = 6, plot = plot_AhK)
    
    
    
    # LEfSe for ANKA_lo vs ANKA_hi
    # ------------------------------
    
    # LEfSE plot
    plot_AlAh <- my_lefserPlot(res_se_AlAh_group, colors = group_colors, 
                              meta_levels=levels(meta_AlAh_group$parasite_strain)) +
      ggtitle("ANKA_lo vs ANKA_hi") + 
      theme(plot.title = element_text(size = 20, face = "bold"))
    plot_AlAh
    
    # Save plot 
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_ANKA_lo_vs_ANKA_hi.png"), width = 10, height = 6, plot = plot_AlAh)
    ggsave(paste("./Output_LEfSe/LEfSe_", paste(taxo_rank, mouse, Days, sep="_"), 
                 "_ANKA_lo_vs_ANKA_hi.svg"), width = 10, height = 6, plot = plot_AlAh)
    

    
  } # End of else if mouse strain
  
  
} # End of function




