

# --------------------------------------------------
# Title : my_lefserPlot
# --------------------------------------------------
# Date : 2021, Jul 02
# Adaptor : Jean-Christophe LONE
# Type : Function
# --------------------------------------------------
# Objective : Plot LEfSe
# Modified from the lefser package
# --------------------------------------------------
# Input(s) : 
# df : data frame lefse results 
# colors : character (color of the groups)
# meta_levels : character (levels of the groups to be associated with colors)
# --------------------------------------------------
# Output(s) : 
# plots : ggplot2 object (plot can be modified)
# --------------------------------------------------



my_lefserPlot <- function (df, colors, meta_levels){
  
  # Debug 
  # trim.names <- FALSE
  # df <- res_se_NK_group
  # colors <-  c("forestgreen", "red")
  # meta_levels <- levels(meta_NK_group$parasite_strain)

  
  # Add groups and levels to 0 and 1 in the initla function
  Group <- ifelse(df$scores > 0, 
                  paste(meta_levels[2]), # One
                  paste(meta_levels[1])) # Zero
  df$Group <- as.factor(Group)


  # Plot the results
  plt <- ggplot(df, aes(reorder(Names, scores), scores)) + 
    ylab("LDA SCORE (log 10)") + 
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_text(size = 20, face = "bold"), 
          axis.text.y = element_text(vjust = 0.7, size = 20, face = "bold"), 
          axis.text.x = element_text(vjust = 0.7, size = 20, face = "bold"), 
          plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text = element_text(color = "black", face = "bold", size = 20),
          legend.title = element_text(color = "black", face = "bold", size = 20)) + 
    geom_bar(stat = "identity", aes(fill = Group), color= "black") + 
    geom_hline(yintercept = 0, colour="black", linetype = "longdash") +
    scale_fill_manual(values = colors) +
    ylim(-4, 4) +
    coord_flip()

  # Get the plot out!
  return(plt)
  
} # End of function





