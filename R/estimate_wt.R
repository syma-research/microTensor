# estimate phi parameter for the quasi-likelihood and 
# corresponding weight tensor for each observation
estimate_wt <- function(X, Yhat) {
  
  Xsum_m <- apply(X, c(2, 3), sum)
  
  expY <- exp(Yhat)
  expYsum_m <- apply(expY, c(2, 3), sum)
  phat <- vapply(seq_len(dim(X)[1]),
                 function(i) expY[i, , ] / expYsum_m,
                 matrix(0.0, nrow = dim(X)[2], ncol = dim(X)[3]))
  phat <- aperm(phat, perm = c(3, 1, 2))
  
  Xhat <- vapply(seq_len(dim(X)[1]),
                 function(i) phat[i, , ] * Xsum_m,
                 matrix(0.0, nrow = dim(X)[2], ncol = dim(X)[3]))
  Xhat <- aperm(Xhat, perm = c(3, 1, 2))
  var_obs <- (X - Xhat)^2
  var_exp <- Xhat * (1 - phat)
  
  ratio <- apply(var_obs, c(2, 3), sum) / apply(var_exp, c(2, 3), sum)
  phi <- mean((ratio - 1) / (Xsum_m - 1), na.rm = TRUE)
  
  wt <- 1 + (Xsum_m - 1) * phi
  return(list(wt = wt,
              phi = phi))
}
