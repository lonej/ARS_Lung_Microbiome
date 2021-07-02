

# --------------------------------------------------
# Title : Adonis_Loop
# --------------------------------------------------
# Date : 2021, Jul 02
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
# PcoA plots : PCoA representation
# stats_adonis : Statistic Table
# vect : Principal Coordinates
# --------------------------------------------------



Adonis_Loop <- function(data, culture, Days, 
                        mouse, parasite, taxo_rank, distance, graph){
  
  
  # Debugging tools
  # rm(list=ls())
  # mouse <- "DBA"
  # culture <- 'IMM'
  # Days <- 5
  # graph <- "Yes"
  # taxo_rank <- "OTU"
  # parasite <- "NAlAhS"
  # distance <- "bray"
  # data <- read.csv2(paste("./Rarefied_Dataset/Rare_exp_",taxo_rank,"_Table.csv",
  #                                  sep=""), check.names=FALSE)
  # data <- data %>% dplyr::select(-Microbiome)
  # data <- data %>%
  #   dplyr::filter( !(mouse_cage %in%
  #                   c("NAl", "NAh", "NS", "AlS", "AhS", "AlAh", "NAlAhS")))

  
  

  
  
  
  # --------------------------------------------------
  # Choose Data subset for PERMANOVA
  # --------------------------------------------------
  
  # Add row names
  row.names(data) <- data$ID
  
  # Select genotype
  sub_Table <- data %>% filter(mouse_strain==mouse)

  
  
  # Day 5 or 7 or both (357)
  # --------------------------------------------------
  if (Days == 357) {
    sub_Table <- sub_Table %>% filter(day_postInfection != 12, culture_type==culture)
    
    # Otherwise choose Mouse strain, culture type and ONE DAY
  } else {
    sub_Table <- sub_Table %>% filter(day_postInfection == Days, culture_type==culture)
  }
  
  
  # Then choose Parasite comparison
  # --------------------------------------------------

    # NI, ANKA_lo, ANKA_hi and SMAC
  if (parasite == "NAlAhS") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | 
                                        parasite_strain=="ANKA_lo" | 
                                        parasite_strain=="ANKA_hi" | 
                                        parasite_strain=="SMAC")
    comparison <- "Global"
    
    # NI vs ANKA low
  } else if (parasite == "NAl") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | parasite_strain=="ANKA_lo")
    comparison <- "NI_vs_ANKA_lo"
    
  } else if (parasite == "NAh") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="NI" | parasite_strain=="ANKA_hi")
    comparison <- "NI_vs_ANKA_hi"
    
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
    
    # ANKA low vs SMAC
  } else if (parasite == "AlS") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="ANKA_lo" | parasite_strain=="SMAC")
    comparison <- "ANKA_lo_vs_SMAC"
    
    # ANKA vs SMAC
  } else if (parasite == "AhS") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="ANKA_hi" | parasite_strain=="SMAC")
    comparison <- "ANKA_hi_vs_SMAC" 
    
    # ANKA vs SMAC
  } else if (parasite == "AlAh") {
    sub_Table <- sub_Table %>% filter(parasite_strain=="ANKA_lo" | parasite_strain=="ANKA_hi")
    comparison <- "ANKA_lo_vs_ANKA_hi" 
    
  } # End of if parasite == NAS 
  
  
  
  
  # --------------------------------------------------
  # PERMANOVA using adonis (vegan)
  # --------------------------------------------------
  
  # Interger Table for PERMANOVA
  sub_Table_int <- sub_Table %>% dplyr::select(-c(ID, parasite_strain, mouse_strain, 
                                                  day_postInfection, culture_type, 
                                                  mouse_cage))
  
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
  
  
  # --------------------------------------------------
  # Principal Coordinate Analysis (PCoA)
  # --------------------------------------------------
  
  # Only build PCoA plots when all 3 parasites are being compared
  if (graph=="Yes" & (parasite == "NAlAhS")){
    
    
    # Eigenvalues evaluation
    # --------------------------------------------------
    
    vect <- as.data.frame(res_betadisp$vectors)
    eig <- res_betadisp$eig
    eig[(eig<0)] <- 0 # Negative eigenvalues are constraint to 0
    eigenvalues <- round(eig/sum(eig)*100, 1)
    
    
    # DBA mouse
    # --------------------------------------------------
    
    if (mouse == "DBA") {
      Groups <- sub_Table %>% dplyr::select(parasite_strain) %>% dplyr::pull() %>% as.factor()
      Groups <- factor(Groups, levels = c("NI", "ANKA_lo", "ANKA_hi", "SMAC"))
      group_colors <- c("black", "orange", "firebrick3", "grey70")
      dot_shape <- c(22, 15, 15, 15)

    
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
  
    } # End if mouse == DBA
    
  } # End if (graph=="Yes" & (parasite == "NAlAhS"))
  
  return(list(stats_adonis=stats_adonis,
              vect=res_betadisp$vectors))
  
} # End function







