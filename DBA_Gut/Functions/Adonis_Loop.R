

# --------------------------------------------------
# Title : Adonis_Loop
# --------------------------------------------------
# Date : 2021, Apr 07
# Author : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Pairwise Adonis and PCoA graphs
# --------------------------------------------------
# Input(s) : 
# data : OTU table with metadata
# culture : in case of bacterial culture
# Days : 3, 5 or 357 (for day 5 and 7 combined)
# mouse : "C57BL6J" or "DBA
# parasite : character ("NI", "ANKA", "K173", "SMAC")
# taxo_rank : "domain", "phylum", "class", "order", "family", "genus", "OTU"
# distance : character ("bray", "jaccard")
# graph : Boolean ("Yes" "No") to display PCoA graphs 
# --------------------------------------------------
# Output(s) : 
# Pco1 plots : PCoA representation
# OTU_*_Pval.csv : Table of paiwise p-values comparing each parsaite infection
# --------------------------------------------------



Adonis_Loop <- function(data, culture, Days, 
                        mouse, parasite, taxo_rank, distance, graph){
  
  
  # Debugging tools
  # rm(list=ls())
  # mouse <- "C57BL6J"
  # culture <- 'IMM'
  # Days <- 5
  # graph <- "Yes"
  # taxo_rank <- "family"
  # parasite <- "NA"
  # distance <- "bray"
  # graph <- "Yes"
  # data <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)

  
  
  

  
  
  
  # --------------------------------------------------
  # Choose Data subset for PERMANOVA
  # --------------------------------------------------
  
  # Add row names
  row.names(data) <- data$ID
  
  # C57BL6J or DBA or All genotype
  # --------------------------------------------------
  if (mouse=="All") {
    sub_Table <- data
  } else {
    sub_Table <- data %>% dplyr::filter(mouse_strain==mouse)
  }
  
  
  # Day 5 or 7 or both (357)
  # --------------------------------------------------
  if (Days == 357) {
    sub_Table <- sub_Table %>% dplyr::filter(day_postInfection != 12, culture_type==culture)
    
    # Otherwise choose Mouse strain, culture type and ONE DAY
  } else {
    sub_Table <- sub_Table %>% dplyr::filter(day_postInfection == Days, culture_type==culture)
  }
  
  
  # Then choose Parasite comparison
  # --------------------------------------------------
    # NI vs ANKA vs K173
  if (parasite == "NAK") {
    sub_Table <- sub_Table
    comparison <- "All_Three"
    
    # NI vs ANKA vs SMAC
  } else if (parasite == "NAS") {
    sub_Table <- sub_Table
    comparison <- "All_Three"
    
    # NI vs ANKA
  } else if (parasite == "NA") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | parasite_strain=="ANKA")
    comparison <- "NI_vs_ANKA"
    
    # ANKA vs K173
  } else if (parasite == "AK") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="ANKA" | parasite_strain=="K173")
    comparison <- "ANKA_vs_K173"
    
    # NI vs K173
  } else if (parasite == "NK") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | parasite_strain=="K173")
    comparison <- "NI_vs_K173"
    
    # NI vs SMAC
  } else if (parasite == "NS") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | parasite_strain=="SMAC")
    comparison <- "NI_vs_SMAC"
    
    # ANKA vs SMAC
  } else if (parasite == "AS") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="ANKA" | parasite_strain=="SMAC")
    comparison <- "ANKA_vs_SMAC" 
  } # End of if parasite == NAS 
  
  
  
  
  # --------------------------------------------------
  # PERMANOVA using adonis (vegan)
  # --------------------------------------------------
  
  # Interger Table for PERMANOVA
  sub_Table_int <- sub_Table %>% dplyr::select(-c(ID, parasite_strain, mouse_strain, 
                                                  day_postInfection, culture_type))
  
  # Permutation constraints
  permutations <- how(nperm = 999)
  setBlocks(permutations) <- with(sub_Table, sub_Table$day_postInfection)
  
  # PERMANOVA using Bray Curtis 
  res_adonis <- adonis(sub_Table_int ~ parasite_strain,
                       data = sub_Table, 
                       permutations = permutations,
                       method=distance)
  
  
  
  # --------------------------------------------------
  #  Multivariate homogeneity of groups dispersions  (Betadisper)
  # --------------------------------------------------
  
  # Metadata Table
  sub_Table_Info <-  sub_Table %>% unite(Info, parasite_strain, mouse_strain,
                                         day_postInfection, culture_type, sep = "_")
  
  # Test of variance equality
  dis_bray <- vegdist(sub_Table_int, method= distance)
  res_betadisp <- betadisper(dis_bray, sub_Table_Info$Info)
  res_anova_beta <- anova(res_betadisp)
  
  
  
  
  
  # --------------------------------------------------
  # Statistics summary
  # --------------------------------------------------
  
  # Data frame of PERMANOVA statistics
  stats_adonis <- cbind.data.frame(Mouse_strain=mouse, 
                                   Culture=culture, 
                                   Days = Days,
                                   Compare = comparison,
                                   Pval = res_adonis$aov.tab["parasite_strain", "Pr(>F)"],
                                   R2 = res_adonis$aov.tab["parasite_strain", "R2"],
                                   Betadisp = res_anova_beta$`Pr(>F)`[1])
  
  # stats_adonis
  
  
  
  
  
  # --------------------------------------------------
  # Principal Coordinate Analysis (PCoA)
  # --------------------------------------------------
  
  # Only build PCoA plots when all 3 parasites are being compared
  if (graph=="Yes" & (parasite == "NAS" | parasite == "NAK")){
    
    
    # Eigenvalues evaluation
    # --------------------------------------------------
    
    vect <- as.data.frame(res_betadisp$vectors)
    eig <- res_betadisp$eig
    eig[(eig<0)] <- 0 # Negative eigenvalues are constraint to 0
    eigenvalues <- round(eig/sum(eig)*100, 1)
    
    
    # What is the genotype considered?
    # --------------------------------------------------

      # Only C57BL6J
    if (mouse == "C57BL6J") {
      Groups <- sub_Table %>% select(parasite_strain) %>% dplyr::pull() %>% as.factor()
      Groups <- factor(Groups, levels = c("NI", "ANKA", "K173"))
      group_colors <- c("black", "snow3", "firebrick1")
      dot_shape <- c(21, 19, 19)
      if (Days==357){
        Groups <- sub_Table %>% unite(PD, parasite_strain, day_postInfection) %>% dplyr::pull(PD) %>% as.factor()
        Groups <- factor(Groups, levels = c("NI_5", "ANKA_5", "K173_5",
                                            "NI_7", "ANKA_7", "K173_7"))
        group_colors <- rep(c("black", "snow3", "firebrick1"),2)
        dot_shape <- c(21, 19, 19, 17, 17, 17)
      }
      
      # Only DBA
    } else if (mouse == "DBA") {
      Groups <- sub_Table %>% select(parasite_strain) %>% dplyr::pull() %>% as.factor()
      Groups <- factor(Groups, levels = c("NI", "ANKA", "SMAC"))
      group_colors <- c("black", "firebrick3", "grey70")
      dot_shape <- c(22, 15, 15)
      if (Days==357){
        Groups <- sub_Table %>% unite(PD, parasite_strain, day_postInfection) %>% dplyr::pull(PD) %>% as.factor()
        Groups <- factor(Groups, levels = c("NI_5", "ANKA_5", "SMAC_5",
                                            "NI_7", "ANKA_7", "SMAC_7"))
        group_colors <- rep(c("black", "firebrick2", "grey70"),2)
        dot_shape <- c(22, 15, 15, 17, 17, 17)
      }
    }
    

    
    
    
    # Finally the PCoA plot
    # --------------------------------------------------
    
    p_bray <- ggplot(data = vect, aes(x = PCoA1, y = PCoA2, group=Groups)) +
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
    p_bray
    
    # Save the rebellion, save the dream
    ggsave(paste("./Output_Adonis/", 
                 paste(taxo_rank, mouse, culture, parasite, Days, distance, sep="_"),
                 "_Axis1-2.png", sep=""), width = 10, height = 6, plot = p_bray)
    ggsave(paste("./Output_Adonis/",
                 paste(taxo_rank, mouse, culture, parasite, Days, distance, sep="_"),
                 "_Axis1-2.svg", sep=""), width = 10, height = 6, plot = p_bray)
  
    
  }
  
  return(stats_adonis)
  
}







