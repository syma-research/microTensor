# main function used to fit microTensor
microTensor <- function(X, R,
                        weighted = TRUE,
                        ortho_m = TRUE,
                        nn_t = TRUE,
                        control = list()) {
  # check X
  
  control <- do.call(control_microTensor, control)
  if(!is.null(control$debug_dir)) {
    dir.create(control$debug_dir, showWarnings = TRUE, recursive = TRUE)
    control$debug_dir <- normalizePath(control$debug_dir)
  }
  
  Yhat <- array(0, dim(X))
  m_prev <- matrix(ncol = 0,
                   nrow = dim(X)[1])
  l_fit_unweighted <- list()
  wt <- array(1, dim(X)[c(2, 3)])
  
  # Initialize
  if(control$verbose)
    message("Initializing...")
  if(control$init == "cp") {
    fit_init <- cp_wrapper(X = X, R = R, nn_t = nn_t)
  } else if (control$init == "ctf") {
    fit_init <- ctf(X = X, R = R, control = list(verbose = FALSE))
  } else if (control$init == "zero") {
    fit_init <- list(lambda = rep(1, R),
                     m = matrix(1 / sqrt(dim(X)[1]), 
                                nrow = dim(X)[1], 
                                ncol = R),
                     s = matrix(1 / sqrt(dim(X)[2]), 
                                nrow = dim(X)[2], 
                                ncol = R),
                     t = matrix(1 / sqrt(dim(X)[3]), 
                                nrow = dim(X)[3], 
                                ncol = R))
  } else {
    stop("Unrecognized initial value!")
  }
  
  # First fit without weighting
  if(control$verbose)
    message("Performing unweighted fitting...")
  for(r in seq(1, R)) {
    if(control$verbose)
      message("Fitting for R = ", r, "...")
    i_fit <- solve_oneR(X = X,
                        Yhat = Yhat,
                        wt = wt,
                        m_prev = m_prev,
                        init = list(lambda = fit_init$lambda[r],
                                    m = fit_init$m[, r, drop = FALSE],
                                    s = fit_init$s[, r, drop = FALSE],
                                    t = fit_init$t[, r, drop = FALSE]),
                        ortho_m = ortho_m,
                        nn_t = nn_t,
                        gamma = control$gamma, L_init = control$L_init,
                        abs_tol = control$abs_tol, rel_tol = control$rel_tol,
                        maxit = control$maxit,
                        debug_dir = paste_debugDir(control$debug_dir,
                                                   "/unweighted_R_", r), 
                        verbose = control$verbose) 
    Yhat <- Yhat + get_Yhat(i_fit)
    m_prev <- cbind(m_prev, matrix(i_fit$m, ncol = 1))
    l_fit_unweighted[[r]] <- i_fit
  }
  
  # Weighted fitting
  if(weighted) {
    if(control$verbose)
      message("Performing weighted fitting...")
    
    Yhat_unweighted <- Yhat
    Yhat <- array(0, dim(X))
    m_prev <- matrix(ncol = 0,
                     nrow = dim(X)[1])
    l_fit_weighted <- list()
    wt_estimate <- estimate_wt(X = X, Yhat = Yhat)
    wt <- wt_estimate$wt
    
    for(r in seq(1, R)) {
      if(control$verbose)
        message("Fitting for R = ", r, "...")
      
      i_fit <- solve_oneR(X = X,
                          Yhat = Yhat,
                          wt = wt,
                          m_prev = m_prev,
                          init = l_fit_unweighted[[r]],
                          ortho_m = ortho_m,
                          nn_t = nn_t,
                          gamma = control$gamma, L_init = control$L_init,
                          abs_tol = control$abs_tol, rel_tol = control$rel_tol,
                          maxit = control$maxit,
                          debug_dir = paste_debugDir(control$debug_dir,
                                                     "/weighted_R_", r), 
                          verbose = control$verbose) 
      Yhat <- Yhat + get_Yhat(i_fit)
      m_prev <- cbind(m_prev, matrix(i_fit$m, ncol = 1))
      l_fit_weighted[[r]] <- i_fit
    }
  }
  
  # reorder lambda
  results_unweighted <- extract_fits(l_fit_unweighted)
  if(!weighted)
    return(results_unweighted)
  
  results_weighted <- extract_fits(l_fit_weighted)
  return(c(results_weighted,
           list(wt_estimate = wt_estimate,
                unweighted = results_unweighted)))
}

# helper for creating control parameters
control_microTensor <- function(
  gamma = 2, 
  L_init = 2,
  abs_tol = 1e-4,
  rel_tol = 1e-4,
  maxit = 1000,
  init = "cp",
  debug_dir = NULL,
  verbose = TRUE
) {
  list(gamma = gamma,
       L_init = L_init,
       abs_tol = abs_tol,
       rel_tol = rel_tol,
       maxit = maxit,
       init = init,
       debug_dir = debug_dir,
       verbose = verbose)
}
