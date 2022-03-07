#' Tensor decomposition for longitudinal microbiome data
#' 
#' \code{microTensor} takes as input a three-dimensional array of longitudinal
#' microbiome sequence count observations (\code{X}), and performs dimension
#' reduction by approximating the observations with \code{R} factors, each
#' characterized by the outer product of a rank-1 microbe-, subject-, and 
#' time-specific loadings.
#' 
#' \code{control} should be provided as a named list, with (some of) the 
#' following components:
#' \describe{
#' \item{gamma}{gamma parameter for the inner-loop gradient descent algorithm.}
#' \item{L_init}{initial L parameter for the inner-loop gradient descent 
#' algorithm.}
#' \item{abs_tol}{absolute tolerance for algorithm convergence.}
#' \item{rel_tol}{relative tolerance for algorithm convergence.}
#' \item{maxit}{maximum allowed number of outer algorithm iterations.}
#' \item{init}{initialization strategy. \code{"cp"} (default) is the canonical 
#' polyadic decomposition. \code{"ctf"} is initialization with the 
#' compositional tensor factorization method proposed in Martino et al, Nature
#' Biotechnology, 2021. \code{"zero"} is initializing at all zero values for the
#' parameters (all samples are uniform).}
#' \item{phi_iter}{non-negative integer indicating if additional iterations over
#' the dispersion parameter \code{phi} should be performed for weighted 
#' quasi-likelihood estimation. The default is \code{0} which indicates no
#' additional iterations. A positive integer indicates the maximum additional
#' iterations.}
#' \item{debug_dir}{character string indicating the directory to store 
#' intermediate fitted results, for debugging purposes.}
#' \item{verbose}{indicates if verbose fitting progress information should be
#' printed. Default to \code{TRUE}.}
#' }
#' 
#' @param X three-dimensional array of longitudinal microbial sequencing counts.
#' First dimension corresponds to microbes, second subjects, and third time points.
#' @param R the total number of factors (dimensions) to decompose \code{X} into.
#' @param weighted indicates if quasi-likelihood weighted estimates should be
#' performed. This is strongly encouraged, and is default to \code{TRUE}.
#' @param ortho_m indicates if the identified microbial loadings should be
#' constrained to be orthonormal. Orthonormal microbial loadings provide better
#' interpretations for variance decomposition among the identified factors. 
#' Default to \code{TRUE}.
#' @param nn_t indicates if the time loadings should be constrained to be 
#' non-negative. Non-negative time loadings can provide more interpretable 
#' effects on the longitudinal trend effects in the identified factors. Default
#' to \code{TRUE}.
#' @param control a named list of additional control parameters. See details.
#'
#' @return a list, with the following components:
#' \describe{
#' \item{results}{
#' list of fitted parameters. If \code{weighted} is set to \code{TRUE}, this 
#' returns the quasi-likelihood weighted estimates, otherwise the unweighted
#' estimates are returned. The components include:
#' \describe{
#' \item{lambda}{vector of length \code{R} for the fitted singular value for
#' each factor.
#' }
#' \item{m}{matrix of \code{R} columns for the fitted microbial loadings.
#' }
#' \item{s}{matrix of \code{R} columns for the fitted subject loadings.
#' }
#' \item{t}{matrix of \code{R} columns for the fitted time point loadings.
#' }
#' \item{fitting}{list of additional model fitting information, including the
#' fitted objective function value, number of outer iterations, and algorithm
#' convergence information.
#' }
#' }
#' }
#' \item{wt_estimate}{
#' list of fitted quasi-likelihood re-weighting parameters, including
#' \describe{
#' \item{wt}{three-dimensional array for the per-sample weights. This has the
#' same dimensionality as \code{X}.
#' }
#' \item{phi}{estimated model dispersion parameter.
#' }
#' }
#' If \code{weighted} is \code{FALSE}, then this component is empty.
#' }
#' }

#' @export
#' @author Siyuan Ma, \email{syma.research@gmail.com}
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
    
    wt_estimate <- estimate_wt(X = X, Yhat = Yhat)
    wt <- wt_estimate$wt
    Yhat <- array(0, dim(X))
    m_prev <- matrix(ncol = 0,
                     nrow = dim(X)[1])
    l_fit_weighted <- list()
    
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
    
    if(control$phi_iter) {
      i_outer_iter <- 1
      wt_estimate_new <- 
        wt_estimate <- 
        estimate_wt(X = X, Yhat = Yhat)
      wt <- wt_estimate$wt
      Yhat <- array(0, dim(X))
      m_prev <- matrix(ncol = 0,
                       nrow = dim(X)[1])
      
      while(i_outer_iter <= control$phi_iter) {
        if(control$verbose) 
          message("Fitting additional phi iteration ", i_outer_iter, "...")
        
        l_fit_weighted_new <- list()
        for(r in seq(1, R)) {
          if(control$verbose)
            message("Fitting for R = ", r, "...")
          i_fit <- 
            solve_oneR(
              X = X,
              Yhat = Yhat,
              wt = wt,
              m_prev = m_prev,
              init = l_fit_weighted[[r]],
              ortho_m = ortho_m,
              nn_t = nn_t,
              gamma = control$gamma, L_init = control$L_init,
              abs_tol = control$abs_tol, rel_tol = control$rel_tol,
              maxit = control$maxit,
              debug_dir = 
                paste_debugDir(
                  control$debug_dir,
                  "/weighted_iter_", i_outer_iter,
                  "_R_", r), 
              verbose = control$verbose)
          
          Yhat <- Yhat + get_Yhat(i_fit)
          m_prev <- cbind(m_prev, matrix(i_fit$m, ncol = 1))
          l_fit_weighted_new[[r]] <- i_fit
        }
        
        # check convergence
        wt_estimate_new <- 
          estimate_wt(X = X, Yhat = Yhat)
        if(abs(wt_estimate_new$phi - wt_estimate$phi) < control$abs_tol | 
           abs(wt_estimate_new$phi - wt_estimate$phi) / 
           (abs(wt_estimate$phi) + control$abs_tol) < control$rel_tol) {
          break
        }
        
        wt <- wt_estimate$wt
        Yhat <- array(0, dim(X))
        m_prev <- matrix(ncol = 0,
                         nrow = dim(X)[1])
        l_fit_weighted <- l_fit_weighted_new
        wt_estimate <- wt_estimate_new
        i_outer_iter <- i_outer_iter + 1
      }
      
      l_fit_weighted <- l_fit_weighted_new
    }
  }
  
  # reorder lambda
  if(!weighted)
    return(list(results = extract_fits(l_fit_unweighted),
                wt_estimate = NULL))
  
  return(list(results = extract_fits(l_fit_weighted),
              wt_estimate = wt_estimate))
}

control_microTensor <- function(
  gamma = 2, 
  L_init = 2,
  abs_tol = 1e-4,
  rel_tol = 1e-4,
  maxit = 1000,
  init = "cp",
  phi_iter = 0,
  debug_dir = NULL,
  verbose = TRUE
) {
  list(gamma = gamma,
       L_init = L_init,
       abs_tol = abs_tol,
       rel_tol = rel_tol,
       maxit = maxit,
       init = init,
       phi_iter = phi_iter,
       debug_dir = debug_dir,
       verbose = verbose)
}
