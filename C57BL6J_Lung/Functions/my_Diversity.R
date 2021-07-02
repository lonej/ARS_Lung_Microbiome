
# --------------------------------------------------
# Title : my_Diversity
# --------------------------------------------------
# Date : 2021, Jul 02
# Author : Jean-Christophe LONE
# Type : Function
# Version : v3
# --------------------------------------------------
# Objective : Compute Shannon and Chao1 diversity
# --------------------------------------------------
# Input(s) : 
# table : OTU table with metadata
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
# mouse : "C57BL6J" or "DBA
# Days : 3, 5 or 357 (for day 5 and 7 combined)
# --------------------------------------------------
# Output(s) : 
# Shannon diveristy plot : graphics
# Chao1 diversity plots : graphics
# --------------------------------------------------





my_Diversity <- function(table, taxo_rank, Days, mouse){
  

  # Clear memory
  # rm(list=ls())
  # taxo_rank <- "OTU"
  # Days <- 5
  # mouse <- "C57BL6J"
  # table <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # table <- table %>% select(-Microbiome)

  
  # Handle Days and Genotype
  # --------------------------------------------------
  
  # Mouse condition
  if (mouse=="C57BL6J") {
    sub_Table <- table %>% dplyr::filter(mouse_strain==mouse)
  } 
  
  # Days condition
  if (Days==5) {
    sub_Table <- sub_Table %>% filter(day_postInfection==Days)
  } else if (Days==7) {
    sub_Table <- sub_Table %>% filter(day_postInfection==Days)
  } 

  
  # Shannon Diversity
  # --------------------------------------------------
  
  # Add ID as row names and Merge Parasite and mouse strain infos
  sub_Table <- sub_Table %>% column_to_rownames("ID") 
  
  # Disjunction between informations and integers
  sub_Table_int <- sub_Table %>% 
    select(-c(parasite_strain, day_postInfection, culture_type, mouse_strain))
  sub_Table_info <- sub_Table %>% 
    select(c(parasite_strain, day_postInfection, culture_type, mouse_strain))
  
  # Shannon index
  shannon <- diversity(sub_Table_int, index = "shannon", MARGIN = 1, base = exp(1)) 
  shannon <- cbind.data.frame(sub_Table_info, shannon)
  

  # Stats
  KW_shannon <- kruskal.test(shannon ~ parasite_strain, data = shannon)
  Wil_shannon <- pairwise.wilcox.test(shannon$shannon, shannon$parasite_strain,
                       p.adjust.method = "holm")
  tmp <- Wil_shannon$p.value %>% melt() %>% filter(value !="NA")
  colnames(tmp) <- c("group1", "group2", "Pval")
  tmp <- cbind.data.frame(tmp)
  
  
  
  # Graphical parameters for ggplot
  # --------------------------------------------------
  
  if (mouse == "C57BL6J") {
    shannon <- shannon %>% dplyr::rename(Groups=parasite_strain)
    shannon$Groups <- factor(shannon$Groups, levels = c("NI", "K173", "ANKA"))
    my_comparison <- list( c("NI", "K173"),  c("K173", "ANKA"), c("ANKA", "NI"))
    group_colors <- c("white",  "firebrick1", "snow3")
  } 
  
  
  # Plot Shannon
  # --------------------------------------------------


  
  # Plot Shannon diversity
  shannon_p <- ggplot(shannon, aes(x=Groups, y=shannon, fill=Groups)) + 
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center') +
    scale_fill_manual(values=group_colors) +
    theme(axis.text.x = element_text(angle = 0, hjust=0.5, size=25,
                                     colour = "black", face = "plain"),
          axis.text.y = element_text(size=15,  colour = "black")) +
    ylab("Alpha Diversity (Shannon)") + 
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title = element_text(color = "black", face = "plain", size = 25),
      legend.text = element_text(color = "black", face = "plain", size = 25),
      axis.title.y = element_text(size = 25, face = "plain", margin = unit(c(0,6,0,0), "mm")),
      axis.title.x = element_blank(),
      axis.text=element_text(size=25)
    )
  shannon_p
  
  # Save the Shannon Graph
  ggsave(paste("./Output_Diversity/Shannon_",  mouse, "_", taxo_rank, "_D", Days,".png", sep=""),
         width = 8, height = 6, plot = shannon_p)
  
  ggsave(paste("./Output_Diversity/Shannon_",  mouse, "_", taxo_rank, "_D", Days,".svg", sep=""),
         width = 8, height = 6, plot = shannon_p)
  
  
  
  
  # Chao Diversity
  # --------------------------------------------------
  
  # Function to run Chao diversity on all columns
  Chao_Richness <- function(Dataset){
    # Initialisation
     results <- data.frame(matrix(ncol = 1, nrow = 2))
    # Compute H for each column
    for (i in seq(ncol(Dataset))) {
      c1 <- chao1(Dataset[,i], taxa.row = TRUE)
      c2 <- chao2(Dataset[,i], taxa.row = TRUE)
      # Get the results in same data frame
      iter <- rbind.data.frame(c1, c2)
      colnames(iter) <- colnames(Dataset)[i]
      results <- cbind.data.frame(results, iter)
      row.names(results) <- c("Chao1", "Chao2")}
    # Output
    results <- results[-c(2),-c(1)]
    return(results)
  }
  
  
  # Chao index computation
  Chao <- Chao_Richness(t(sub_Table_int)) %>% t
  Chao <- cbind.data.frame(sub_Table_info, Chao)
  
  # Stats
  KW_Chao <- kruskal.test(Chao1 ~ parasite_strain, data = Chao)
  Wil_Chao <- pairwise.wilcox.test(Chao$Chao, Chao$parasite_strain,
                                      p.adjust.method = "BH")

  
  if (mouse == "C57BL6J") {
    Chao <- Chao %>% dplyr::rename(Groups=parasite_strain)
    Chao$Groups <- factor(Chao$Groups, levels = c("NI", "K173", "ANKA"))
    my_comparison <- list( c("NI", "K173"),  c("K173", "ANKA"), c("ANKA", "NI"))
    group_colors <- c("white",  "firebrick1", "snow3")
  } 
  
  
  
  # Plot Shannon diversity
  # --------------------------------------------------
  
  # Shannon ggplot
  chao_p <-  ggplot(Chao, aes(x=Groups, y=Chao1, fill=Groups)) + 
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center') +
    scale_fill_manual(values=group_colors) +
    theme(axis.text.x = element_text(angle = 0, hjust=0.5, size=25,
                                     colour = "black", face = "plain"),
          axis.text.y = element_text(size=15,  colour = "black")) +
    ylab("Alpha Diversity (Chao)") + 
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title = element_text(color = "black", face = "plain", size = 25),
      legend.text = element_text(color = "black", face = "plain", size = 25),
      axis.title.y = element_text(size = 25, face = "plain", margin = unit(c(0,6,0,0), "mm")),
      axis.title.x = element_blank(),
      axis.text=element_text(size=25)
    )
  chao_p
  
  # Save the Chao
  ggsave(paste("./Output_Diversity/Chao_",  mouse, "_", taxo_rank, "_D", Days,".png", sep=""),
         width = 8, height = 6, plot = chao_p)
  ggsave(paste("./Output_Diversity/Chao_",  mouse, "_", taxo_rank, "_D", Days,".svg", sep=""),
         width = 8, height = 6, plot = chao_p)
  
  
  # Save the diveristy scores
  if (mouse=="C57BL6J"){
    diver <- rbind.data.frame(
      cbind.data.frame(Method="Shannon", Compare="All", Pval=KW_shannon$p.value),
      cbind.data.frame(Method="Shannon", Compare="NI vs ANKA", Pval=Wil_shannon$p.value[2,1]),
      cbind.data.frame(Method="Shannon", Compare="NI vs K173", Pval=Wil_shannon$p.value[2,2]),
      cbind.data.frame(Method="Shannon", Compare="ANKA vs K173", Pval=Wil_shannon$p.value[1,1]),
      cbind.data.frame(Method="Chao1", Compare="All", Pval=KW_Chao$p.value),
      cbind.data.frame(Method="Chao1", Compare="NI vs ANKA", Pval=Wil_Chao$p.value[2,1]),
      cbind.data.frame(Method="Chao1", Compare="NI vs K173", Pval=Wil_Chao$p.value[2,2]),
      cbind.data.frame(Method="Chao1", Compare="ANKA vs K173", Pval=Wil_Chao$p.value[1,1])
    )
  } 

  write.csv2(diver,paste("./Output_Diversity/Stats_",  mouse, "_", taxo_rank, "_D", Days,".csv", sep="") )
  
}
  
  
  







