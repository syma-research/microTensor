# wrapper for the canonical polyadic decomposition, used
# for microTensor
cp_wrapper <- function(X,
                       Yhat = array(0, dim(X)),
                       R = 1,
                       centering_m = TRUE,
                       m_proj = NULL,
                       nn_t = TRUE) {
  # if there are missing values in X
  X_norm_fill <- X
  X_norm_fill[is.na(X_norm_fill)] <- 0
  fill_val <- min(setdiff(X_norm_fill, 0))  / 2
  X_norm_fill <- apply(X_norm_fill + fill_val, c(2, 3),
                       function(x) x / sum(x))
  
  # suppress progress bar
  message_cp <- capture.output(
    fit_cp <- rTensor::cp(
      tnsr = rTensor::as.tensor(log(X_norm_fill) - Yhat), 
      num_components = R))
  
  lambda <- fit_cp$lambdas
  m <- fit_cp$U[[1]]
  s <- fit_cp$U[[2]]
  t <- fit_cp$U[[3]]
  
  # transform m
  if(centering_m)
    m <- apply(m, 2, function(x) x - mean(x))
  if(!is.null(m_proj)) {
    if(R != 1)
      stop("m_proj is only valid for when R = 1!")
    if(centering_m)
      stop("m centering should be disabled when m_proj is provided!")
    m <- m_proj %*% m
  }
  
  for(r in seq(1, R)) {
    # sign flip
    if(lambda[r] < 0) {
      m[, r] <- -m[, r]
      lambda[r] <- -lambda[r]
    }
    if(s[order(-abs(s[, r]))[1], r] < 0) {
      m[, r] <- -m[, r]
      s[, r] <- -s[, r]
    }
    if(t[order(-abs(t[, r]))[1], r] < 0) {
      m[, r] <- -m[, r]
      t[, r] <- -t[, r]
    }
    if(nn_t)
      t[t[, r] < 0, r] <- 0
    
    # normalize
    norms <- sqrt(c(sum(m[, r]^2),
                    sum(s[, r]^2),
                    sum(t[, r]^2)))
    m[, r] <- m[, r] / norms[1]
    s[, r] <- s[, r] / norms[2]
    t[, r] <- t[, r] / norms[3]
    lambda[r] <- lambda[r] * prod(norms)
  }
  
  lambda_order <- order(-lambda)
  lambda <- lambda[lambda_order]
  m <- m[, lambda_order, drop = FALSE]
  s <- s[, lambda_order, drop = FALSE]
  t <- t[, lambda_order, drop = FALSE]
  
  return(list(lambda = lambda,
              m = m,
              s = s, 
              t = t,
              fitting = list(
                convergence = 1 - fit_cp$conv
              )))
}

# wrapper for the canonical polyadic decomposition, used
# for CTF
cp_wrapper_ctf <- function(X,
                           R = 1) {
  # if there are missing values in X
  X_fill <- X
  X_fill[is.na(X_fill)] <- 0
  
  # suppress progress bar
  message_cp <- capture.output(
    fit_cp <- rTensor::cp(
      tnsr = rTensor::as.tensor(X_fill), 
      num_components = R))
  
  lambda <- fit_cp$lambdas
  m <- fit_cp$U[[1]]
  s <- fit_cp$U[[2]]
  t <- fit_cp$U[[3]]
  
  for(r in seq(1, R)) {
    # sign flip
    if(lambda[r] < 0) {
      m[, r] <- -m[, r]
      lambda[r] <- -lambda[r]
    }
    if(s[order(-abs(s[, r]))[1], r] < 0) {
      m[, r] <- -m[, r]
      s[, r] <- -s[, r]
    }
    if(t[order(-abs(t[, r]))[1], r] < 0) {
      m[, r] <- -m[, r]
      t[, r] <- -t[, r]
    }
    
    # normalize
    norms <- sqrt(c(sum(m[, r]^2),
                    sum(s[, r]^2),
                    sum(t[, r]^2)))
    m[, r] <- m[, r] / norms[1]
    s[, r] <- s[, r] / norms[2]
    t[, r] <- t[, r] / norms[3]
    lambda[r] <- lambda[r] * prod(norms)
  }
  
  lambda_order <- order(-lambda)
  lambda <- lambda[lambda_order]
  m <- m[, lambda_order, drop = FALSE]
  s <- s[, lambda_order, drop = FALSE]
  t <- t[, lambda_order, drop = FALSE]
  
  return(list(lambda = lambda,
              m = m,
              s = s, 
              t = t,
              fitting = list(
                convergence = 1 - fit_cp$conv
              )))
}
