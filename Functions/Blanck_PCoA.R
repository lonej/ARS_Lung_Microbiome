
# --------------------------------------------------
# Title : Adonis_Loop
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Show Blanck have different "microbiome"
# --------------------------------------------------
# Input(s) : 
# data : OTU table with metadata
# mouse2rm : List of cages ti remove
# --------------------------------------------------
# Output(s) : 
# DBA_Blanck_PCoA.png : PCoA Plot with Blank
# --------------------------------------------------




Blanck_PCoA <- function(data, mouse2rm){
  
  # Debug
  # rm(list=ls())
  # taxonomy_rank <- "OTU"
  # mouse2rm <- c("D7_2", "D7_7" , "D7_10")
  # exp_OTU_Table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_WithBlanck_",taxonomy_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # data <-  exp_OTU_Table %>% dplyr::select(-Microbiome)
  
  
  
  # Remove outlier cages
  data <- data %>% filter( !(mouse_cage %in% mouse2rm)) %>% 
    column_to_rownames("ID")
  
  
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
  
  
  
  # Plot one parasite strain
  
  Groups <- data %>% dplyr::select(parasite_strain) %>% dplyr::pull() %>% as.factor()
  Groups <- factor(Groups, levels = c("NI", "ANKA", "SMAC", "Sentinel", "DBA_Blanck"))
  group_colors <- c("black", "firebrick3", "grey70", "forestgreen", "blue")
  dot_shape <- c(22, 15, 15, 15, 15)
  
  
  plot1 <- ggplot(data = vect, aes(x = PCoA1, y = PCoA2, group=Groups)) +
    geom_point(aes(color=Groups, shape=Groups), size=5) +
    scale_shape_manual(values=dot_shape) +
    scale_color_manual(values=group_colors) +
    # geom_text_repel(aes(label=row.names(vect), color=Groups), size = 1.75, max.overlaps = Inf) +
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
  
  
  ggsave(paste("./Checkpoints_Sentinels_and_Blanck/DBA_Blanck_PCoA.png", sep=""),
         width = 10, height = 6, plot = plot1)
  
  return(plot1)
}




