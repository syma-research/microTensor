# wrapper for performing PCA
pca <- function(X, R,
                control = list()) {
  dims <- dim(X)
  X <- t(apply(X, 1, as.vector))
  
  fill_val <- min(setdiff(X, 0), na.rm = TRUE)  / 2
  X_norm_fill <- apply(X + fill_val, 2,
                       function(x) x / sum(x))
  ind_missing <- apply(X, 2, function(x) all(is.na(x)))
  
  fit_prcomp <- prcomp(t(log(X_norm_fill[, !ind_missing])))
  
  m <- fit_prcomp$rotation[, seq(1, R)]
  s_full <- matrix(NA, nrow = ncol(X_norm_fill), ncol = R)
  s_full[!ind_missing, ] <- fit_prcomp$x[, seq(1, R)]
  
  lambda <- apply(s_full, 2, function(x) sqrt(sum(x^2, na.rm = TRUE)))
  s_full <- apply(s_full, 2, function(x) x / sqrt(sum(x^2, na.rm = TRUE)))
  
  s <- apply(s_full, 2, function(r_s) {
    matrix(r_s, nrow = dims[2], ncol = dims[3]) %>% 
      apply(1, mean, na.rm = TRUE)
  })
  
  return(list(lambda = lambda,
              m = m,
              s = s,
              s_full = s_full))
}
