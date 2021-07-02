

# --------------------------------------------------
# Title : Clean_N_Rarefied
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : 
# Clean, Check and Merge Metadata and OTU Table
# --------------------------------------------------
# Input(s) : 
# min_reads : integer (minimal number of read to consider)
# rare_curve : boolean (display a graph or not)
# --------------------------------------------------
# Output(s) :
# Read_count_OTU_Histogram.png : Raster graphical (repartition of reads per sample)
# Samples_over_min_reads.csv : Table (sample over the threshold)
# Samples_under_min_reads.csv : Table (sample below the threshold)
# Root_removal.csv : Table (OTU identified as Root)
# Chloroplast_removal.csv : Table (OTU identified as chloroplast)
# Rarefaction_curve.png : Raster graphics (Rarefaction curve)
# Rare_exp_WithBlanck_ : Rarified Dataset with blank
# Rare_exp_ : Rarified Dataset without blanck
# --------------------------------------------------


Clean_N_Rarefied = function(min_reads, rare_curve){
 
  
  # Debugging
  # rm(list=ls())
  # min_reads <- 1000
  # rare_curve <- FALSE
  
  # ----------------------------------------------
  # Clean Metadata
  # --------------------------------------------------
  
  
  # Download Meta data
  metadata <- read.csv2("./Dataset/Metadata.csv", header=T, sep = ";", na.strings = "")

  # Just rename sample.id to ID
  metadata <- metadata %>% filter(!is.na(sample.id)) %>% dplyr::rename(ID = sample.id)
  
  # Replace NAs in column by their "Empty" name
  metadata <- metadata %>% replace_na(list(parasite_strain = "NI",
                                           mouse_strain = "No_Strain",
                                           day_postInfection = 0,
                                           culture_type = "Lung",
                                           mouse_cage = "No_Cage"))
  
  # dplyr::select the informative variables
  metadata <- metadata %>% dplyr::select(ID, Microbiome, parasite_strain, mouse_strain, 
                                  day_postInfection, culture_type, mouse_cage)
  
  
  # --------------------------------------------------
  # Clean OTU Table
  # --------------------------------------------------
  
  # Download Datasets 
  OTU_Table <- read.delim("./Dataset/OTU_Table.txt", header = TRUE, sep = "\t", dec = ".")
  Taxonomy <- read.delim("./Dataset/Taxonomy.txt", header = TRUE, sep = "\t", dec = ".")
  
  # Rename X.OTU.ID colomn by OTU
  OTU_Table <- OTU_Table %>% dplyr::rename(OTU = X.OTU.ID)
  
  # Remove "_paired_M05485" from colnames
  colnames(OTU_Table) <- OTU_Table %>% colnames() %>% gsub("_paired_M05485", "", .)
  
  # Use OTU names as row names
  OTU_Table <- OTU_Table %>% column_to_rownames("OTU") 


  
  
  # --------------------------------------------------
  # Remove samples based on reads number
  # --------------------------------------------------

  # Number of reads per sample
  png("./Checkpoints_Cleaning/Read_count_OTU_Histogram.png", units="in", width=7, height=5, res=1200)
  OTU_Table %>%  colSums() %>% 
    hist(breaks=50, main="Lung - DBA", 
         xlab="Number of reads",
         ylab = "Number of samples", border="blue", 
         col="dodgerblue", las=2,  xaxt = "n", cex=0.5)  %>% 
    axis(side=1, at=seq(0,2000000, 10000),  labels=seq(0,2000000,10000)) 
  dev.off()
  
  # Remove Low sequences samples
  reads_per_sample <- OTU_Table %>% t %>% rowSums() %>% as.data.frame()
  colnames(reads_per_sample) <- "Count"
  over_min_reads <- reads_per_sample %>% filter(Count>min_reads) 
  under_min_reads <- reads_per_sample %>% filter(Count<min_reads)
  OTU_Table <- OTU_Table %>% dplyr::select(all_of(row.names(over_min_reads)))
  
  # Save csv file 
  write.csv2(over_min_reads, "./Checkpoints_Cleaning/Samples_over_min_reads.csv")
  write.csv2(under_min_reads, "./Checkpoints_Cleaning/Samples_under_min_reads.csv")

  
  
  
  # --------------------------------------------------
  # Remove unclassified_Root and Chloroplast
  # --------------------------------------------------
  
  # Split the "OTU1;size=55875" and keep OTU only
  # The warning sign is because it discarded the second column (the ;size=123)
  Taxonomy <- Taxonomy %>% separate(X, c("OTU", NA))
  
  # Remove unclassified_Root and chloroplast
  Root <- Taxonomy[Taxonomy$domain=="unclassified_Root", c(1, 2)]
  Chloro <- Taxonomy[Taxonomy$class=="Chloroplast", c(1, 4)] 
  otu_remove <- c(Root$OTU, Chloro$OTU)
  OTU_Table <- OTU_Table %>% filter( !(row.names(OTU_Table) %in% otu_remove))
  
  # Save the infos w
  write.csv2(Root, "./Checkpoints_Cleaning/Root_removal.csv")
  write.csv2(Chloro, "./Checkpoints_Cleaning/Chloroplast_Removal.csv")

  

  
  
  # --------------------------------------------------
  # Rarefaction
  # --------------------------------------------------

  # Expected"rarefied number of species
  sRare <- rarefy(OTU_Table, min(rowSums(OTU_Table)), se = FALSE, MARGIN = 1)
  
  OTU_Table <- OTU_Table %>% t.data.frame()
  
  # Rarefaction curves
  if (rare_curve == TRUE){
  source("./Functions/quickRareCurve.R")
  png("./Checkpoints_Cleaning/Rarefaction_curve.png", 
      units="in", width=8, height=5, res=1200)
  curveRare <- quickRareCurve(OTU_Table, step = 20, 
                         sample = min(rowSums(OTU_Table)), 
                         col = "blue", cex = 0.6, label = FALSE)
  dev.off()
  }
  
  # Randomly rarefied community data frame
  rare <- rrarefy(OTU_Table, min(rowSums(OTU_Table)))
  

  
  
  
  
  # --------------------------------------------------
  # Add taxonomy infos to OTU Table
  # --------------------------------------------------

  # Add full taxonomy list to OTU_Table
  rare <- t(rare)
  rare <- cbind.data.frame(OTU=row.names(rare), rare)
  Tax_OTU_Table <- inner_join(Taxonomy, rare, by="OTU")
  
  # Checks
  Tax_OTU_Table[Tax_OTU_Table$genus=="unclassified_Root",c(1, 7)] %>% nrow()
  Tax_OTU_Table[Tax_OTU_Table$genus=="unclassified_Bacteria",c(1, 7)] %>% nrow()
  Tax_OTU_Table[Tax_OTU_Table$genus=="chloroplast",c(1, 7)] %>% nrow()
  anti_join(Taxonomy, rare, by="OTU") %>% nrow()
  anti_join(rare, Taxonomy, by="OTU") %>% nrow()

  
  
  
  # --------------------------------------------------
  # Create a Table for each taxonomy rank
  # --------------------------------------------------

  # taxo_rank <- "family"
  
  for (taxo_rank in c("domain", "phylum", "class",
                      "order", "family", "genus", "OTU")){
    
    # Chose Taxonomy level
    OneTax_Table <- Tax_OTU_Table %>% 
      dplyr::select(all_of(taxo_rank), 8:ncol(Tax_OTU_Table)) %>% 
      dplyr::rename(Taxo = all_of(taxo_rank))
    
    # Sum all read from same taxon
    OneTax_Table <- aggregate(OneTax_Table %>% dplyr::select(-Taxo),
                              by = list(OneTax_Table$Taxo),
                              FUN = "sum")
    
    OneTax_Table <- OneTax_Table %>% column_to_rownames('Group.1') %>% t.data.frame


    
    
    # --------------------------------------------------
    # Merging Data and Metadata
    # --------------------------------------------------
    
    # Add ID column
    OneTax_Table <- cbind.data.frame(ID = row.names(OneTax_Table), OneTax_Table)
    
    # Add Metadata information based on ID
    exp_OTU <- inner_join(metadata, OneTax_Table, by="ID")
    
    # Save Dataset with Blanck in it
    write.csv2(exp_OTU, 
               paste("./Rarefied_Dataset/Rare_exp_WithBlanck_",taxo_rank,
                     "_Table.csv", sep=""),  row.names=F)
    
    # Dataset without Blanck
    blanck_list <- exp_OTU %>% 
      filter(mouse_strain=="B6_Blanck" | mouse_strain=="DBA_Blanck") %>% 
      dplyr::select(ID)
    exp_OTU <- exp_OTU %>% dplyr::filter(!(ID %in% blanck_list$ID))

    # Save Dataset
    write.csv2(exp_OTU, paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,
                              "_Table.csv", sep=""),row.names=F)
    
  } # End of taxo_rank for loop 
   
} # End function










