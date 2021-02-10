# ------------------------------------------------------------------------
# -----------------GroupBootPlot------------------------------------------
# ------------------------------------------------------------------------



# -----------------GroupBootPlot------------------------------------------



GroupBootPlot <- function(BootOut, GroupNames, edges.x, edges.y, 
                          labels = TRUE, transparency = 0.15, point.size = 1.5, 
                          line.size = 1, scales = "fixed",
                          legend.position = "none", GroupNamesCheck = FALSE)
   
{
   
   
   ## Use BootTable to prepare a data.frame for plotting from GroupNetworkBoot output
   boot_dat <- BootTable(BootOut)
   
   
   ## rbind edge data to separate sample and boot.mean for plotting
   boot_dat <- rbind(boot_dat, boot_dat)
   boot_dat$mean.type <- rep(c("Sample", "Bootstrap Mean"), each = nrow(boot_dat) / 2) 
   boot_dat$mean <- boot_dat$sample
   boot_dat[boot_dat$mean.type == "Bootstrap Mean",]$mean <- boot_dat[boot_dat$mean.type == "Bootstrap Mean",]$boot.mean
   
   
   
   ## If GroupNamesCheck == FALSE, print name matching to console
   if(!missing(GroupNames) & GroupNamesCheck == TRUE)  {
      GroupMatch <- cbind.data.frame(data.list.names = sort(names(BootOut$sample$network)), 
                                     GroupNames = GroupNames)
      for(i in 1:nrow(GroupMatch))  {cat("data.list ", as.character(GroupMatch[i, 1]), 
                                         " matches GroupNames ", as.character(GroupMatch[i, 2]), "\n")}
   } 
   
   ## Define GroupNames if not indicated
   if(missing(GroupNames)) {GroupNames <- names(BootOut$sample$network)}
   
   # Create GroupNames variable in boot_dat
   boot_dat$gnames <- NA
   
   for(i in 1:length(unique(boot_dat$g)))  {
      
      boot_dat[boot_dat$g == unique(boot_dat$g)[i], "gnames"] <- GroupNames[i]
      
   }
   
   
   ## Rank edges by their summed size across groups
   
   # Create rank data
   edge_rank <- boot_dat[boot_dat$mean.type == "Sample", c("edges", "mean")] %>%
      group_by(edges) %>%
      summarise_all(sum) %>%
      arrange(mean)
   
   # use rank data to inform factor levels
   boot_dat$edges <- factor(boot_dat$edges, levels = edge_rank$edges)
   
   
   
   ## If edges.x & edges.y are defined, subset to unique edge combinations
   
   if(!missing(edges.x) & !missing(edges.y))      {
      
      # get unique edge combinations
      edges <- c(as.vector(outer(edges.x, edges.y, FUN = "paste", sep = "-")),
                 as.vector(outer(edges.y, edges.x, FUN = "paste", sep = "-")))
      
      # Subset to relevant edges only
      boot_dat <- boot_dat[boot_dat$edges %in% edges,]
      
   }
   
   
   ## Create plot
   boot_plot <- ggplot(boot_dat, aes(y = mean, x = edges, group = .data$mean.type))  +
      geom_point(aes(col = .data$mean.type), size = point.size) +
      geom_line(aes(col = .data$mean.type), size = line.size) +
      geom_ribbon(aes(ymin = .data$ci.lb, ymax = .data$ci.ub), alpha = transparency) +
      coord_flip() +
      scale_color_manual(values = c("darkred", "black")) +
      geom_hline(yintercept = 0, col = "black") +
      facet_grid(. ~ gnames, scales = scales) +
      theme_bw() +
      {if(!labels) theme(axis.text.y = element_blank())} +
      theme(legend.position = legend.position,
            legend.title = element_blank(),
            axis.title = element_blank()) 
   
   # Return plot as output
   return(boot_plot)
}


