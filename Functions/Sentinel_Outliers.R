
# --------------------------------------------------
# Title : Adonis_Loop
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Find outlier cages with sentinels
# --------------------------------------------------
# Input(s) : 
# data : OTU table with metadata
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
#
# --------------------------------------------------
# Output(s) : 
# PcoA plots : PCoA representation
#
# --------------------------------------------------



Sentinel_Outliers <- function(data, taxonomy_rank){
  
  
  # Debug
  # rm(list=ls())
  # taxonomy_rank <- "OTU"
  # exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxonomy_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # data <-  exp_OTU_Table %>% dplyr::select(-Microbiome)
  
  
  
  
  # Select Variables
  data <- data %>% filter(parasite_strain == "Sentinel") %>% 
    column_to_rownames("ID")
  
  
  
  # --------------------------------------------------
  # PERMANOVA using adonis (vegan)
  # --------------------------------------------------
  
  # Interger Table for PERMANOVA
  sub_Table_int <- data %>% dplyr::select(-c(parasite_strain, mouse_strain, 
                                             day_postInfection, culture_type,
                                             mouse_cage))
  
  # Metadata Table
  sub_Table_Info <-  data %>% unite(Info, parasite_strain, mouse_strain,
                                    day_postInfection, culture_type, sep = "_")
  
  # Permutation constraints
  permutations <- how(nperm = 999)
  setBlocks(permutations) <- with(data, data$day_postInfection)
  
  # Test of variance equality
  dis_bray <- vegdist(sub_Table_int, method= "bray")
  res_betadisp <- betadisper(dis_bray, sub_Table_Info$Info)
  
  # Eigenvalues evaluation
  vect <- as.data.frame(res_betadisp$vectors)
  eig <- res_betadisp$eig
  eig[(eig<0)] <- 0 # Negative eigenvalues are constraint to 0
  eigenvalues <- round(eig/sum(eig)*100, 1)
  
  
  
  # --------------------------------------------------
  # PCoA plot of Sentinels
  # --------------------------------------------------
  
  # Graphical options
  group_colors <- c("forestgreen")
  dot_shape <- c(15)
  
  # Plot one parasite strain
  Groups <- data %>% dplyr::select(parasite_strain) %>% dplyr::pull() %>% as.factor()
  
  plot1 <- ggplot(data = vect, aes(x = PCoA1, y = PCoA2, group=Groups)) +
    geom_point(aes(color=Groups, shape=Groups), size=5) +
    scale_shape_manual(values=dot_shape) +
    scale_color_manual(values=group_colors) +
    geom_text_repel(aes(label=row.names(vect), color=Groups), size = 1.75, max.overlaps = Inf) +
    ylab(paste("PCoA Axis 2 [", eigenvalues[2], "%]", sep="")) +
    xlab(paste("PCoA Axis 1 [", eigenvalues[1], "%]", sep="")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title = element_text(color = "black", face = "plain", size = 20),
      legend.text = element_text(color = "black", face = "plain", size = 20),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title.y = element_text(size = 25, face = "plain", margin = unit(c(0,6,0,0), "mm")),
      axis.title.x = element_text(size = 25, face = "plain", margin = unit(c(6,0,0,0), "mm")),
      axis.text=element_text(size=15)
    )
  plot1
  
  ggsave(paste("./Checkpoints_Sentinels_and_Blanck/DBA_Sentinels_PCoA.png", sep=""),
         width = 10, height = 6, plot = plot1)
  
  
  # --------------------------------------------------
  # Outlier cage list
  # --------------------------------------------------
  
  # Looking for outlier Sentinels for dim.1
  sentinels_outliers <- 
    vect %>% dplyr::filter(PCoA1 < quantile(vect$PCoA1, 0.25)-1.5*IQR(vect$PCoA1) |
                             PCoA1 > quantile(vect$PCoA1, 0.75) + 1.5*IQR(vect$PCoA1)) %>% 
    row.names()
  
  
  return(sentinels_outliers)
  
}




