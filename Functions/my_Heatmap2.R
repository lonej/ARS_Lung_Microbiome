


# --------------------------------------------------
# Title : my_Heatmap2
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Heatmaap of taxons
# --------------------------------------------------
# Input(s) : 
# data : OTU table with metadata
# days : 3, 5 
# max_zero_taxa : Maximal number of zero in taxon allowed
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
#
# --------------------------------------------------
# Output(s) : 
# Heatmaps : Percentage Heatmap
# --------------------------------------------------






my_Heatmap2 <- function(data, days, max_zero_taxa, taxo_rank, mouse){
  
  
  
  # Debugging tools
  # rm(list=ls())
  # days <- 5
  # mouse <- "DBA"
  # max_zero_taxa <- 20
  # taxo_rank <- "family"
  # exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,
  #                                  "_Table.csv", sep=""), check.names=FALSE)
  # data <- exp_OTU_Table %>% dplyr::select(-c("Microbiome", "mouse_cage")) %>%
  #   filter(parasite_strain !="Sentinel")
  
  
  
  
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
    dplyr::select(-c("parasite_strain", "mouse_strain",
                     "day_postInfection", "culture_type")) %>% colSums() 
  
  # Remove empty taxa
  zeros <- addition[addition==0] %>% as.data.frame() %>% row.names()
  data <- data %>% dplyr::select(c("parasite_strain", "mouse_strain",
                                   "day_postInfection", "culture_type"),
                                 !all_of(zeros))
  
  
  
  # --------------------------------------------------
  # Heatmap Data calculations
  # --------------------------------------------------
  
  # Annotation and integer
  data_int <- data %>% dplyr::select(-c("parasite_strain", "mouse_strain",
                                        "day_postInfection", "culture_type"))
  
  data_info <- data %>% dplyr::select(c("parasite_strain", "mouse_strain",
                                        "day_postInfection", "culture_type"))
  
  
  # Count the number of zero per Taxon and remove the high number of zero.
  low_zero_taxa <- colSums(data_int==0) %>% 
    as.data.frame() %>% 
    dplyr::rename(nb_zero=".") %>% 
    dplyr::filter(nb_zero<max_zero_taxa)
  low_zero_taxa_names <- low_zero_taxa %>% rownames()
  
  # Calculate percentages in each sample
  data_per <- sweep(data_int, 1, rowSums(data_int), `/`)
  data_per <- log((data_per*100 + 0.01), 2)
  data_per <- cbind.data.frame(data_info, data_per)
  
  # Select all the low zero taxa
  data_per <- data_per %>% 
    dplyr::select(parasite_strain, mouse_strain, day_postInfection, culture_type,
                  all_of(low_zero_taxa_names))
  
  
  
  # --------------------------------------------------
  # Heatmap Parameters
  # --------------------------------------------------
  
  # Graphical displays
  taxon_text_size <- 5
  mice_text_size <- 7
  width <- 8
  height <- 6
  res <- 1200
  
  # Parameters for pheatmap
  annotation <- data_per %>% dplyr::select("parasite_strain") %>%
    dplyr::rename(Parasite=parasite_strain)
  annotation_colors <- list(
    Parasite = c(NI="white", ANKA_lo ="orange", ANKA_hi ="firebrick3", SMAC="grey70"))
  data_per$parasite_strain <- factor(data_per$parasite_strain, level=c("NI", "SMAC", "ANKA_lo", "ANKA_hi"))
  data_per <- data_per %>% arrange(parasite_strain)
  
  # Change into matrix form to run pheatmap
  mtx_data <- data_per %>%
    dplyr::select(-c("parasite_strain",  "mouse_strain",
                     "day_postInfection", "culture_type")) %>%
    t.data.frame
  
  
  # Breaks scale for colors
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(mtx_data, n = 20)
  
  

  
  
  # --------------------------------------------------
  # Heatmap Plots
  # --------------------------------------------------
  
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
                            filename = paste("./Output_Heatmaps/HM2_",
                                             paste(taxo_rank, mouse, days, sep="_"),
                                             ".png", sep="")
  )
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
                           filename = paste("./Output_Heatmaps/HM2_", 
                                            paste(taxo_rank, mouse, days, sep="_"),
                                            ".pdf", sep=""))
  pheatmap_pdf
  
  
}


