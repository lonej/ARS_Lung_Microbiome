
# --------------------------------------------------
# Title : my_Heatmap
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Heatmap of taxon 
# --------------------------------------------------
# Input(s) : 
# data : OTU table with metadata
# days : 3, 5 or 357 (for day 5 and 7 combined)
# threshold_fold : Integer (number of minimal read per sample)
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
# mouse : "C57BL6J" or "DBA
# --------------------------------------------------
# Output(s) :
# heatmap : Raster graphical .png
# heatmap : Vector graphical .pdf
# --------------------------------------------------



my_Heatmap <- function(data, days, threshold_fold, taxo_rank, mouse){
  
  # Debugging tools
  # rm(list=ls())
  # days <- 7
  # mouse <- "All"
  # mouse <- "C57BL6J"
  # threshold_fold <- 10
  # taxo_rank <- "family"
  # exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,
  #                                  "_Table.csv", sep=""), check.names=FALSE)
  # data <- exp_OTU_Table %>% dplyr::select(-Microbiome)

  
  
  
  # --------------------------------------------------
  # Selection of parameters
  # --------------------------------------------------
  
  # Mouse genotype
  if (mouse=="All"){
    data <- data
  } else if (mouse=="C57BL6J" | mouse=="DBA"){
    data <- data %>% filter(mouse_strain==mouse)
  }
  
  # Data cleaning and selection
  data <- column_to_rownames(data, "ID")
  
  # Care od day values
  if (days==357){
    data <- data
  } else if (days==7 | days==5){
    data <- data %>% filter(day_postInfection==days)
  }

  
  
  
  # --------------------------------------------------
  # Dataframe cleaning
  # --------------------------------------------------
  
  # Total sum of each taxon 
  addition <- data %>% 
    select(-c("parasite_strain", "mouse_strain",
              "day_postInfection", "culture_type")) %>% colSums() 
  
  # Remove empty taxa
  zeros <- addition[addition==0] %>% as.data.frame() %>% row.names()
  data <- data %>% select(c("parasite_strain", "mouse_strain",
                                "day_postInfection", "culture_type"),
                              !all_of(zeros))
  
  
  
  # --------------------------------------------------
  # Heatmap Data calculations
  # --------------------------------------------------
  
  # Annotation and integer
  data_int <- data %>% select(-c("parasite_strain", "mouse_strain",
                                     "day_postInfection", "culture_type"))
    
  data_info <- data %>% select(c("parasite_strain", "mouse_strain",
                                      "day_postInfection", "culture_type"))
  
  # Calculate percentages in each sample
  data_per <- sweep(data_int, 1, rowSums(data_int), `/`)
  data_per <- log((data_per*100 + 0.01), 2)
  data_per <- cbind.data.frame(data_info, data_per)
  
  
  
  # --------------------------------------------------
  # Criterias to show a taxon
  # --------------------------------------------------
  
  # Number of read per taxon in each experimental group
  
  f_NI_B6 <- data %>% filter(parasite_strain=="NI", mouse_strain=="C57BL6J") %>%
    select(-c("parasite_strain",  "mouse_strain",
              "day_postInfection", "culture_type")) %>% colSums()
  
  f_ANKA_B6 <- data %>% filter(parasite_strain=="ANKA", mouse_strain=="C57BL6J") %>%
    select(-c("parasite_strain",  "mouse_strain",
              "day_postInfection", "culture_type")) %>% colSums()
  
  f_K173_B6 <- data %>% filter(parasite_strain=="K173", mouse_strain=="C57BL6J") %>%
    select(-c("parasite_strain",  "mouse_strain",
              "day_postInfection", "culture_type")) %>% colSums()
  
  
  f_compare <- cbind.data.frame(f_NI_B6, f_ANKA_B6, f_K173_B6)
  
  # Count the number of samples in each group of each genotype
  thresh_NI_B6 <- data_info %>% filter(parasite_strain=="NI", mouse_strain=="C57BL6J") %>% nrow()
  thresh_ANKA_B6 <-data_info %>% filter(parasite_strain=="ANKA", mouse_strain=="C57BL6J") %>% nrow()
  thresh_K173_B6 <- data_info %>% filter(parasite_strain=="K173", mouse_strain=="C57BL6J") %>% nrow()


  # Selection of Taxons above a threshold
  if (mouse=="C57BL6J"){
    selected <- f_compare %>% 
      dplyr::filter((f_NI_B6 > threshold_fold*thresh_NI_B6) |
                      (f_ANKA_B6 > threshold_fold*thresh_ANKA_B6) |
                      (f_K173_B6 > threshold_fold*thresh_K173_B6)) %>% row.names()
  } 
  

  
  
  # --------------------------------------------------
  # Heatmap
  # --------------------------------------------------
  
  # Reorder Dataset
  data_per %>% row.names()
  
  
  # Graphical parameters
  annotation <- data_per %>% select("parasite_strain") %>%
    dplyr::rename(Parasite=parasite_strain)
  annotation_colors <- list(
    Parasite = c(NI="white", K173="firebrick1",  ANKA="snow3"))
  data_per$parasite_strain <- factor(data_per$parasite_strain, level=c("NI", "K173", "ANKA"))
  data_per <- data_per %>% arrange(parasite_strain)
  
  
  # Change into matrix form to run pheatmap
  mtx_data <- data_per %>%
    select(-c("parasite_strain",  "mouse_strain",
              "day_postInfection", "culture_type")) %>%
    dplyr::select(all_of(selected)) %>% t.data.frame
  
  
  # Breaks scale for colors
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(mtx_data, n = 100)
  

  # Graphical displays
  if (days == 357){
    taxon_text_size <- 5
    mice_text_size <- 7
    width <- 7
    height <- 5
    res <- 1200
  } else {
    taxon_text_size <- 5
    mice_text_size <- 7
    width <- 8
    height <- 6
    res <- 1200
  }


    
  # Plot Heatmap
  pheatmap_data <- pheatmap(mtx_data,  angle_col = 90,
                            cluster_cols = F, cluster_rows =T,
                            show_colnames = F, show_rownames = T,
                            annotation_col = annotation, 
                            annotation_colors = annotation_colors,
                            fontsize_row = taxon_text_size, 
                            fontsize_col =mice_text_size,
                            color = plasma(length(mat_breaks) - 1),
                            breaks = mat_breaks,
                            cellwidth = width, cellheight = height,
                            filename = paste("./Output_Heatmaps/HM_", 
                                             paste(taxo_rank, mouse, days, sep="_"),
                                             ".png", sep=""))
  pheatmap_data
  
  
  # Save in a vector format (pdf because svg is not permitted)
  pheatmap_pdf <- pheatmap(mtx_data,  angle_col = 90,
                            cluster_cols = F, cluster_rows =T,
                            show_colnames = F, show_rownames = T,
                            annotation_col = annotation, 
                            annotation_colors = annotation_colors,
                            fontsize_row = taxon_text_size, 
                            fontsize_col =mice_text_size,
                            color = plasma(length(mat_breaks) - 1),
                            breaks = mat_breaks,
                            cellwidth = width, cellheight = height,
                            filename = paste("./Output_Heatmaps/HM_", 
                                             paste(taxo_rank, mouse, days, sep="_"),
                                             ".pdf", sep=""))
  pheatmap_pdf
  
} # End of function
  
  
  