
# -----------------BootTable-------------------------------

BootTable <- function(
  BootOut  ## Output of GroupNetworkBoot
) {
  
  ## Define Parameters
  nodes <- colnames(BootOut$sample[[1]][[1]])
  Gn <- length(BootOut$sample[[1]])
  Gnames <- names(BootOut$sample[[1]])
  nboots <- length(BootOut$boot)
  
  ## Create data.frame with all edge combinations
  edges <- do.call(expand.grid, rep(list(nodes), 2))
  edges <- mutate_all(edges, as.character)
  edges <- edges[edges$Var1 != edges$Var2,]
  edges$edges <- apply(cbind(edges$Var1, edges$Var2), 1, function(x) paste(sort(x), collapse="-"))
  edges <- edges[!duplicated(edges$edges),]
  
  ## Create empty vectors to fill
  edges$sample <- NA
  edges$boot.mean <- NA
  edges$ci.lb <- NA
  edges$ci.ub <- NA
  edges$boot.zero <- NA
  
  ## Expand edges data.frame Gn times and save group variable
  edges <- do.call("rbind", replicate(Gn, edges, simplify = FALSE))
  edges$g <- rep(Gnames, each = nrow(edges) / Gn)
  
  ## Run loops to save values
  for(i in 1:nrow(edges)) {
    
    ## Sample Mean
    edges[i, "sample"] = BootOut$sample$network[[edges$g[i]]][edges$Var1[i], edges$Var2[i]]
    
    ## Get Bootstrap vector
    boot_vec <- c()
    
    for(j in 1:nboots)   {
      
      boot_vec[j] <- BootOut$boot[[paste("b", j, sep = "")]][[edges$g[i]]][edges$Var1[i], edges$Var2[i]]
      
    }
    
    ## Get bootstrap values
    edges[i, "boot.mean"] = mean(boot_vec)
    edges[i, c("ci.lb", "ci.ub")] = quantile(boot_vec, probs = c(0.025, 0.975))
    edges[i, "boot.zero"] = sum(boot_vec == 0) / nboots
    edges[i, "boot.pos"] = sum(boot_vec > 0) / nboots
    edges[i, "boot.neg"] = sum(boot_vec < 0) / nboots
    
    ## Remove temporary variables
    rm(boot_vec)
    
  }
  
  ## Compute the number of time bootstrap sample was 0
  
  ## Sort edge data by sample value
  edges <- arrange(edges, .data$g, .data$sample)
  
  ## Return output
  return(edges)
  
  
}


