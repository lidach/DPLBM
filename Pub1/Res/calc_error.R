# calculate relative error function
calc_error <- function(true, est, iters){
  list <- list(abs_error = NA, rel_error = NA, MRE = NA, MARE = NA)
  # absolute error
  for (i in 1:iters){
    list$abs_error[i] <- (est[i] - true[i])
  }
  # relative error
  for(i in 1:iters){
    list$rel_error[i] <- (est[i] - true[i]) / true[i]
  }
  # bias
  list$MRE <- median(list$rel_error, na.rm = TRUE)
  # precision
  list$MARE <- median(abs(list$rel_error), na.rm =TRUE)
  
  return(list)
}