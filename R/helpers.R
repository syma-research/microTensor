# helper functions

# extract fitted parameters and generate output list
extract_fits <- function(l_fit) {
  lambda_order <- order(-vapply(l_fit,
                              function(x) x$lambda,
                              0.0))
  l_fit <- l_fit[lambda_order]
  
  results <- list(lambda = vapply(l_fit,
                                  function(x) x$lambda,
                                  0.0),
                  m = vapply(l_fit,
                             function(x) x$m[, 1],
                             rep(0.0, length = length(l_fit[[1]]$m))),
                  s = vapply(l_fit,
                             function(x) x$s[, 1],
                             rep(0.0, length = length(l_fit[[1]]$s))),
                  t = vapply(l_fit,
                             function(x) x$t[, 1],
                             rep(0.0, length = length(l_fit[[1]]$t))),
                  fitting = list(
                    obj = vapply(l_fit,
                                       function(x) x$obj,
                                       0.0),
                    steps = vapply(l_fit,
                                   function(x) x$steps,
                                   0.0),
                    convergence = vapply(l_fit,
                                         function(x) x$convergence,
                                         0.0)
                  ))
  
  return(results)
}

# create the Y tensor from fitted parameters
get_Yhat <- function(result, R = length(result$lambda)) {
  R <- length(result$lambda)
  l_Y <- 
    lapply(seq(1, R),
           function(r) {
             result$lambda[r] * 
               outer(result$m[, r], 
                     outer(result$s[, r],
                           result$t[, r]))
           })
  Yhat <- 
    Reduce("+", l_Y)
  return(Yhat)
}

# RCLR transform as defined in Martino et al.
rclr <- function(x) {
  # if this observation is purely NA
  if(all(is.na(x)))
    return(x)
  
  if(any(is.na(x)))
    stop("x should have either all missing values, or none!")
  if(all(x == 0))
    stop("x cannot all be zeros!")
  if(any(x < 0))
    stop("Cannot be applied to negative values!")
  
  x[x == 0] <- NA
  # x <- x / sum(x, na.rm = TRUE)
  log(x) - mean(log(x), na.rm = TRUE)
}

# utility for debug directories
paste_debugDir <- function(parent_dir, ...) {
  if(is.null(parent_dir))
    return(NULL)
  else
    return(paste0(parent_dir, ...))
}
