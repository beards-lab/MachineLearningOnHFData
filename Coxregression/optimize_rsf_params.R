optimize_rsf_params <- function(data, param_grid, selected_features) {
  
  best_oob <- Inf
  best_params <- list()
  param_results <- data.frame(
    ntree = numeric(), 
    mtry = numeric(), 
    nodesize = numeric(), 
    nodedepth = character(), 
    splitrule = character(), 
    samptype = character(),
    avg_oob = numeric(), 
    stringsAsFactors = FALSE
  )
  
  # Ensure selected_features includes time and event columns
  selected_features <- c(selected_features, "time", "event")
  
  # Subset the data to include only the selected features
  data <- data[, selected_features]
  
  total_combinations <- length(param_grid$ntree) * length(param_grid$mtry) * 
    length(param_grid$nodesize) * length(param_grid$nodedepth) * 
    length(param_grid$splitrule) * length(param_grid$samptype)
  progress <- 0
  
  # Iterate over all combinations in the parameter grid
  for (ntree in param_grid$ntree) {
    for (mtry in param_grid$mtry) {
      for (nodesize in param_grid$nodesize) {
        for (nodedepth in param_grid$nodedepth) {
          for (splitrule in param_grid$splitrule) {
            for (samptype in param_grid$samptype) {
              
              oob_scores <- numeric()
              
              for (i in 1:25) {
                # Handle "NoLimit" case for nodedepth explicitly
                if (nodedepth == "NoLimit") {
                  rf_model <- rfsrc(Surv(time, event) ~ ., data = data, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                    splitrule = splitrule, samptype = samptype)
                } else {
                  rf_model <- rfsrc(Surv(time, event) ~ ., data = data, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                    nodedepth = as.numeric(nodedepth), splitrule = splitrule, samptype = samptype)
                }
                
                # Record OOB error
                oob_scores <- c(oob_scores, rf_model$err.rate[length(rf_model$err.rate)])
              }
              
              avg_oob <- mean(oob_scores)
              param_results <- rbind(param_results, data.frame(
                ntree = ntree, mtry = mtry, nodesize = nodesize, 
                nodedepth = nodedepth, splitrule = splitrule, 
                samptype = samptype, avg_oob = avg_oob,
                stringsAsFactors = FALSE
              ))
              
              if (avg_oob < best_oob) {  # Select best parameters based on OOB error
                best_oob <- avg_oob
                best_params <- list(ntree = ntree, mtry = mtry, nodesize = nodesize, nodedepth = nodedepth, splitrule = splitrule, samptype = samptype)
              }
              
              # Update progress
              progress <- progress + 1
              cat("Progress: ", round(100 * progress / total_combinations, 2), "% - Avg OOB: ", avg_oob, "\n")
              
              # Plot OOB error evolution
              plot(param_results$avg_oob, type = "l", col = "blue", xlab = "Parameter Set", ylab = "Average OOB Error",
                   main = "OOB Error Progress")
              abline(h = best_oob, col = "red", lty = 2)
            }
          }
        }
      }
    }
  }
  
  # Final OOB error plot
  ggplot(param_results, aes(x = 1:nrow(param_results), y = avg_oob)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = best_oob, color = "red", linetype = "dashed") +
    labs(title = "OOB Error Across Parameter Combinations", x = "Parameter Set", y = "Average OOB Error")
  
  return(list(best_params = best_params, param_results = param_results))
}