# main function used to fit microTensor
microTensor_iter <- function(X, R,
                             ortho_m = TRUE,
                             nn_t = TRUE,
                             max_outer_iter = 100,
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
  
  ll_fit <- list()
  i_outer_iter <- 1
  while(i_outer_iter <= max_outer_iter) {
    # fit microTensor
    l_fit <- list()
    for(r in seq(1, R)) {
      if(control$verbose)
        message("Fitting for R = ", r, "...")
      if(i_outer_iter == 1) {
        r_fit <- solve_oneR(X = X,
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
      } else {
        r_fit <- solve_oneR(X = X,
                            Yhat = Yhat,
                            wt = wt,
                            m_prev = m_prev,
                            init = ll_fit[[i_outer_iter - 1]]$l_fit[[r]],
                            ortho_m = ortho_m,
                            nn_t = nn_t,
                            gamma = control$gamma, L_init = control$L_init,
                            abs_tol = control$abs_tol, rel_tol = control$rel_tol,
                            maxit = control$maxit,
                            debug_dir = paste_debugDir(control$debug_dir,
                                                       "/unweighted_R_", r), 
                            verbose = control$verbose) 
      }
      
      Yhat <- Yhat + get_Yhat(r_fit)
      m_prev <- cbind(m_prev, matrix(r_fit$m, ncol = 1))
      l_fit[[r]] <- r_fit
    }
    
    # estimate weighting
    wt_estimate <- estimate_wt(X = X, Yhat = Yhat)
    wt <- wt_estimate$wt
    ll_fit[[i_outer_iter]] <- list(l_fit = l_fit,
                                   phi = wt_estimate$phi,
                                   Yhat = Yhat)
    
    # re-initialize for the next loop
    Yhat <- array(0, dim(X))
    m_prev <- matrix(ncol = 0,
                     nrow = dim(X)[1])
    
    # check convergence
    if(i_outer_iter >=2) {
      if(abs(ll_fit[[i_outer_iter]]$phi - ll_fit[[i_outer_iter - 1]]$phi) < 
         control$abs_tol |
         abs(ll_fit[[i_outer_iter]]$phi - ll_fit[[i_outer_iter - 1]]$phi) / 
         (ll_fit[[i_outer_iter - 1]]$phi + control$abs_tol) < control$rel_tol)
        break
    }
    i_outer_iter <- i_outer_iter + 1
  }
  
  return(ll_fit)
}