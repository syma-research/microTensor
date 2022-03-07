ctf <- function(X, R,
                control = list()) {
  control <- do.call(control_microTensor, control)
 
  X <- apply(X, c(2, 3), rclr)
  Yhat <- array(0, dim(X))
  l_fit <- list()
  
  for(r in seq(1, R)) {
    if(control$verbose)
      message("Fitting for R = ", r, "...")
    i_fit <- solve_oneR_ctf(X = X - Yhat,
                            abs_tol = control$abs_tol, rel_tol = control$rel_tol,
                            maxit = control$maxit,
                            debug_dir = paste_debugDir(control$debug_dir,
                                                       "/R_", r), 
                            verbose = control$verbose) 
    Yhat <- Yhat + get_Yhat(i_fit)
    l_fit[[r]] <- i_fit
  }
  
  # reorder lambda
  results <- extract_fits(l_fit)
  return(results)
}
