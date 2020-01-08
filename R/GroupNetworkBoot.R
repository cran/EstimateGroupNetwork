# ------------------------------------------------------------------------
# -----------------Bootstrap Functions for Joint Network------------------
# ------------------------------------------------------------------------



# -----------------GroupNetworkBoot---------------------------------------

GroupNetworkBoot <- function(
   data_list, 
   groupNetwork, 
   nboots = 100, 
   bootSeed,
   ...)     {
   
   
   #----------------Random seed for replicability --------
   if(!missing(bootSeed)) set.seed(bootSeed)
   
   
   # ---------------Check simplifyOutput argument--------
   
   if(is.null(groupNetwork$network[[1]])) {
      stop("Please set simplifyOutput=FALSE in EstimateGroupNetwork of your input network.\nOtherwise this function will not run.")
   }
   
   
   # ---------------Save important variables-------------
   G = length(data_list) # number of samples/groups
   Gnames = names(data_list) # Get sample data.frame names
   Gn = data.frame(G = 1:G) # Data.frame including sample sizes
   Gn$n = sapply(data_list, nrow) 
   
   nvar = ncol(groupNetwork$network[[1]])
   edges = colnames(groupNetwork$network[[1]])
   
   tracker = round(seq(from = 0, to = 100, length.out = nboots), 1)
   
   
   start_time = Sys.time()
   
   
   # --------------- Check data -------------------------
   if(any(sapply(data_list, ncol) != nvar))
   {
      stop("Error: all datasets (data_list) must include the same variables of the original network (groupNetwork)")
   }

   # ---------------Define empty arguments---------------

   labels = colnames(groupNetwork$network[[1]])

   # ---------------Create bootstrapping list------------
   output <- list(data = data_list,
                  sample = groupNetwork,
                  boot = list())
   
   
   # ---------------Run bootstrapping--------------------
   for(i in 1:nboots)      {
      
      # Bootstrap samples
      boot_dat = list()
      
      for(j in 1:G)     {
         
         ## Save group
         boot_dat[[j]] <- data_list[[j]][sample(1:Gn$n[j], size = Gn$n[j], replace = TRUE),]
         
      }
      
      ## Name boot_dat data.frames
      names(boot_dat) <- Gnames
      
      # ---------------Run EstimateGroupNetwork-------
 
      boot_network <- 
         EstimateGroupNetwork(boot_dat,
                              inputType = "list.of.dataframes", 
                              simplifyOutput = FALSE,
                              ...
                              )
   
      
      ## Save bootstrap network edges to output list
      output$boot[[length(output$boot) + 1]] <- boot_network$network
      names(output$boot) <- paste("b", 1:i, sep = "")
      
      
      ## Print tracker to show progress
      
      # Calculate remaining time
      cur_time <- Sys.time()
      time_completed <- as.numeric(cur_time - start_time)
      time_remaining <- round((time_completed / i) * (nboots - i))
      hours_remaining <- round(time_remaining / 60)
      minutes_remaining <- time_remaining - hours_remaining * 60
      
      # Print tracker
      cat(paste("\r", tracker[i], "% (~", hours_remaining, " Hours ", 
                minutes_remaining, " Minutes remaining)        ", sep = ""))
      
      
      # Remove temporary variables
      rm(boot_dat)
      rm(boot_network)
      rm(time_remaining)
      rm(time_completed)
      rm(hours_remaining)
      rm(minutes_remaining)
      rm(cur_time)
      
   }
   
   # Return Output
   return(output)
   
}

