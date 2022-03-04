# fitting function for each R

solve_oneR <- function(X, 
                       Yhat = array(0, dim(X)), 
                       wt = array(1, dim(X)[c(2, 3)]),
                       m_prev = matrix(ncol = 0,
                                       nrow = dim(X)[1]),
                       init = NULL,
                       ortho_m = TRUE,
                       nn_t = TRUE,
                       gamma = 1.25, L_init = 1,
                       abs_tol = 1e-3,
                       rel_tol = 1e-3,
                       maxit = 50,
                       debug_dir = NULL,
                       verbose = TRUE) {
  
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(!all(dim(X)[c(2, 3)] == dim(wt)))
    stop("Dimensions of X and wt do not agree.")
  if(nrow(m_prev) != dim(X)[1])
    stop("Dimensions of X and m_prev do not agree.")
  if(!is.null(debug_dir)) {
    dir.create(debug_dir, showWarnings = TRUE, recursive = TRUE)
    debug_dir <- normalizePath(debug_dir)
  }
  
  # check missingness pattern of X 
  ## FIXME
  
  # constraint for m vector based on specification
  m_constraint <- matrix(rep(sqrt(1/dim(X)[1]),
                             dim(X)[1]),
                         ncol = 1)
  if(ortho_m)
    m_constraint <- cbind(m_constraint, m_prev)
  if(max(abs(t(m_constraint) %*% m_constraint) - 
         diag(rep(1, ncol(m_constraint)))) > 1e-14)
    stop("Projection matrix for m is not orthonormal.")
  m_proj <- diag(rep(1, dim(X)[1])) -
    m_constraint %*% t(m_constraint)
  
  if(is.null(init)) {
    if(verbose)
      message("Initializing...")
    
    # Initialize using Canonical Polyadic decomposition
    fit_cp <- cp_wrapper(X = X, Yhat = Yhat, R = 1, 
                         centering = FALSE,
                         m_proj = m_proj,
                         nn_t = nn_t)
    lambda_init <- fit_cp$lambda
    vm_init <- fit_cp$m[, 1]
    vs_init <- fit_cp$s[, 1]
    vt_init <- fit_cp$t[, 1]
  } else {
    lambda_init <- init$lambda
    vm_init <- init$m[, 1]
    vs_init <- init$s[, 1]
    vt_init <- init$t[, 1]
    
    # make sure that lambda_init is not too large
    max_logY <- max(lambda_init * outer(vm_init, outer(vs_init, vt_init)))
    if(max_logY > 500)
      lambda_init <- lambda_init * 500 / max_logY
  }
  
  f_init <- negLogLik_cpp(X = X,
                      Yhat = Yhat,
                      wt = wt,
                      lambda = lambda_init,
                      vm = vm_init,
                      vs = vs_init,
                      vt = vt_init,
                      normalize = TRUE)
  
  vm_old <- list(vm = vm_init,
                 L = L_init,
                 f = f_init)
  vs_old <- list(vs = vs_init,
                 L = L_init,
                 f = f_init)
  vt_old <- list(vt = vt_init,
                 L = L_init,
                 f = f_init)
  lambda_old <- list(lambda = lambda_init,
                     L = L_init,
                     f = f_init)
  i_step <- 1
  if(!is.null(debug_dir)) {
    l_oneFit <- list(vm = vm_old,
                     vs = vs_old,
                     vt = vt_old,
                     lambda = lambda_old)
    save(l_oneFit, file = paste0(debug_dir, "/fit_", i_step, ".RData"))
  }
  
  # optimize with gradient descent
  if(verbose)
    message("Fitting gradient descent...")
  while(TRUE) {
    if(i_step %% 10 == 0)
      if(verbose)
        message("Iteration ", i_step)
    # update m
    vm_new <- update_m(vm_old = vm_old$vm, X = X, Yhat = Yhat, wt = wt, 
                       m_proj = m_proj,
                       lambda = lambda_old$lambda, vs = vs_old$vs, 
                       vt = vt_old$vt,
                       gamma = gamma, L_old = vm_old$L)
    
    # update s
    vs_new <- update_s(vs_old = vs_old$vs, X = X, Yhat = Yhat, wt = wt, 
                       lambda = lambda_old$lambda, vm = vm_new$vm, 
                       vt = vt_old$vt,
                       gamma = gamma, L_old = vs_old$L)
    
    # update t
    vt_new <- update_t(vt_old = vt_old$vt, X = X, Yhat = Yhat, wt = wt, 
                       nn_t = nn_t,
                       lambda = lambda_old$lambda, vm = vm_new$vm, 
                       vs = vs_new$vs,
                       gamma = gamma, L_old = vt_old$L)
    
    # update lambda
    lambda_new <- update_lambda(lambda_old = lambda_old$lambda, X = X, 
                                Yhat = Yhat, wt = wt,
                                vm = vm_new$vm, vs = vs_new$vs, vt = vt_new$vt,
                                gamma = gamma, L_old = lambda_old$L)
    if(lambda_new$f > vt_new$f)
      stop("Something went wrong!") ##FIXME?
    
    # check for convergence
    if(abs(lambda_old$f - lambda_new$f) < abs_tol &
       abs(lambda_old$f - lambda_new$f) / (abs(lambda_old$f) + abs_tol) < rel_tol) {
      convergence <- 0
      break
    }
    
    vm_old <- vm_new
    vs_old <- vs_new
    vt_old <- vt_new
    lambda_old <- lambda_new
    i_step <- i_step + 1
    
    if(!is.null(debug_dir)) {
      l_oneFit <- list(vm = vm_old,
                       vs = vs_old,
                       vt = vt_old,
                       lambda = lambda_old)
      save(l_oneFit, file = paste0(debug_dir, "/fit_", i_step, ".RData"))
    }
    
    if(i_step > maxit) {
      convergence <- 1
      break
    }
  }
  
  # identifiability check
  if(lambda_new$lambda < 0) {
    vm_new$vm <- -vm_new$vm
    lambda_new$lambda <- -lambda_new$lambda
  }
  if(vs_new$vs[order(-abs(vs_new$vs))[1]] < 0) {
    vm_new$vm <- -vm_new$vm
    vs_new$vs <- -vs_new$vs
  }
  # renormalization
  norms <- sqrt(c(sum(vm_new$vm^2),
                  sum(vs_new$vs^2),
                  sum(vt_new$vt^2)))
  vm_new$vm <- vm_new$vm / norms[1]
  vs_new$vs <- vs_new$vs / norms[2]
  if(norms[3] != 0) { ## FIXME for when vt get shrinked to all zeros!
    vt_new$vt <- vt_new$vt / norms[3]
    lambda_new$lambda <- lambda_new$lambda * prod(norms)
  } else {
    lambda_new$lambda <- lambda_new$lambda * prod(norms[c(1, 2)])
  }
  
  return(list(m = matrix(vm_new$vm, ncol = 1),
              s = matrix(vs_new$vs, ncol = 1),
              t = matrix(vt_new$vt, ncol = 1),
              lambda = lambda_new$lambda,
              obj = lambda_new$f,
              steps = i_step,
              convergence = convergence))
}

solve_oneR_ctf <- function(X,
                           abs_tol = 1e-3,
                           rel_tol = 1e-3,
                           maxit = 50,
                           debug_dir = NULL,
                           verbose = TRUE) {
  
  if(!is.null(debug_dir)) {
    dir.create(debug_dir, showWarnings = TRUE, recursive = TRUE)
    debug_dir <- normalizePath(debug_dir)
  }
  
  
  if(verbose)
    message("Initializing...")
  
  # Initialize using Canonical Polyadic Decomposition
  fit_cp <- cp_wrapper_ctf(X = X, R = 1)
  lambda <- fit_cp$lambda
  vm_init <- fit_cp$m[, 1]
  vs_init <- fit_cp$s[, 1]
  vt_init <- fit_cp$t[, 1]
  
  f_init <- mean((X - lambda * outer(vm_init, outer(vs_init, vt_init)))^2,
                 na.rm = TRUE)
  
  vm_old <- list(vm = vm_init,
                 f = f_init)
  vs_old <- list(vs = vs_init,
                 f = f_init)
  vt_old <- list(vt = vt_init,
                 f = f_init)
  i_step <- 1
  if(!is.null(debug_dir)) {
    l_oneFit <- list(vm = vm_old,
                     vs = vs_old,
                     vt = vt_old)
    save(l_oneFit, file = paste0(debug_dir, "/fit_", i_step, ".RData"))
  }
  
  # optimize with gradient descent
  if(verbose)
    message("Fitting gradient descent...")
  while(TRUE) {
    if(i_step %% 10 == 0)
      if(verbose)
        message("Iteration ", i_step)
    # update m
    covariates_m <- lambda * outer(vs_old$vs, vt_old$vt)
    new_vm <- vapply(seq(1, dim(X)[1]),
                     function(i) {
                       if(!all(is.na(X[i, , ]))) {
                         if(sum((covariates_m * !is.na(X[i, , ]))^2) == 0)
                           return(0)
                         else
                           return(sum(covariates_m * X[i, , ], na.rm = TRUE) /
                                    sum((covariates_m * !is.na(X[i, , ]))^2)) 
                       } else
                         return(0)
                     },
                     0.0)
    vm_new <- list(vm = new_vm,
                   f = mean((X - lambda * outer(new_vm, outer(vs_old$vs, vt_old$vt)))^2,
                            na.rm = TRUE))
    
    # update s
    covariates_s <- lambda * outer(vm_new$vm, vt_old$vt)
    new_vs <- vapply(seq(1, dim(X)[2]),
                     function(j) {
                       sum(covariates_s * X[, j, ], na.rm = TRUE) /
                         sum((covariates_s * !is.na(X[, j, ]))^2)
                     },
                     0.0)
    vs_new <- list(vs = new_vs,
                   f = mean((X - lambda * outer(vm_new$vm, outer(new_vs, vt_old$vt)))^2,
                            na.rm = TRUE))
    
    # update t
    covariates_t <- lambda * outer(vm_new$vm, vs_new$vs)
    new_vt <- vapply(seq(1, dim(X)[3]),
                     function(k) {
                       sum(covariates_t * X[, , k], na.rm = TRUE) /
                         sum((covariates_t * !is.na(X[, , k]))^2)
                     },
                     0.0)
    vt_new <- list(vt = new_vt,
                   f = mean((X - lambda * outer(vm_new$vm, outer(vs_new$vs, new_vt)))^2,
                            na.rm = TRUE))
    
    # check for convergence
    if(abs(vt_new$f - vt_old$f) < abs_tol &
       abs(vt_new$f - vt_old$f) / (abs(vt_old$f) + abs_tol) < rel_tol) {
      convergence <- 0
      break
    }
    
    vm_old <- vm_new
    vs_old <- vs_new
    vt_old <- vt_new
    i_step <- i_step + 1
    
    if(!is.null(debug_dir)) {
      l_oneFit <- list(vm = vm_old,
                       vs = vs_old,
                       vt = vt_old)
      save(l_oneFit, file = paste0(debug_dir, "/fit_", i_step, ".RData"))
    }
    
    if(i_step > maxit) {
      convergence <- 1
      break
    }
  }
  
  # identifiability check
  if(vs_new$vs[order(-abs(vs_new$vs))[1]] < 0) {
    vm_new$vm <- -vm_new$vm
    vs_new$vs <- -vs_new$vs
  }
  # renormalization
  norms <- sqrt(c(sum(vm_new$vm^2),
                  sum(vs_new$vs^2),
                  sum(vt_new$vt^2)))
  vm_new$vm <- vm_new$vm / norms[1]
  vs_new$vs <- vs_new$vs / norms[2]
  vt_new$vt <- vt_new$vt / norms[3]
  lambda <- lambda * prod(norms)
  
  return(list(m = matrix(vm_new$vm, ncol = 1),
              s = matrix(vs_new$vs, ncol = 1),
              t = matrix(vt_new$vt, ncol = 1),
              lambda = lambda,
              obj = vt_new$f,
              steps = i_step,
              convergence = convergence))
}
