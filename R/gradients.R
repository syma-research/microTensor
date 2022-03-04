# These are the gradients for each parameter in the 
# constrained coordinate gradient descent algorithm

gradient_m <- function(X, Yhat, wt, lambda, vm, vs, vt,
                       normalize = TRUE) {
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(dim(X)[1] != length(vm))
    stop("Dimensions of X and vm do not agree.")
  if(dim(X)[2] != length(vs))
    stop("Dimensions of X and vs do not agree.")
  if(dim(X)[3] != length(vt))
    stop("Dimensions of X and vt do not agree.")
  
  st <- outer(vs, vt)
  Xsum_m <- apply(X, c(2, 3), sum)
  expY <- exp(Yhat + lambda * outer(vm, outer(vs, vt)))
  expYsum_m <- apply(expY, c(2, 3), sum)
  term_prod <- Xsum_m / expYsum_m
  
  val <- - vapply(seq(1, length(vm)),
                  function(i) {
                    sum((lambda * st *
                           (X[i, , ] - term_prod * expY[i, , ])) /
                          wt, 
                        na.rm = TRUE)
                  },
                  0.0)
  if(normalize)
    val <- val / 
    median(Xsum_m, na.rm = TRUE) / # normalize by median library size
    sum(!is.na(Xsum_m)) * # normalize by sample size
    median(wt, na.rm = TRUE) # normalize by median weight
  
  return(val)
}

gradient_s <- function(X, Yhat, wt, lambda, vm, vs, vt,
                       normalize = TRUE) {
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(dim(X)[1] != length(vm))
    stop("Dimensions of X and vm do not agree.")
  if(dim(X)[2] != length(vs))
    stop("Dimensions of X and vs do not agree.")
  if(dim(X)[3] != length(vt))
    stop("Dimensions of X and vt do not agree.")
  
  mt <- outer(vm, vt)
  Xsum_m <- apply(X, c(2, 3), sum)
  expY <- exp(Yhat + lambda * outer(vm, outer(vs, vt)))
  expYsum_m <- apply(expY, c(2, 3), sum)
  
  val <- - vapply(seq(1, length(vs)),
                  function(j) {
                    sum((apply(lambda * mt * X[, j, ], 2, sum) -
                           Xsum_m[j, ] *
                           apply(expY[, j, ] * lambda * mt, 2, sum) /
                           expYsum_m[j, ]) / 
                          wt[j, ],
                        na.rm = TRUE)
                  },
                  0.0)
  if(normalize)
    val <- val / 
    median(Xsum_m, na.rm = TRUE) / # normalize by median library size
    sum(!is.na(Xsum_m)) * # normalize by sample size
    median(wt, na.rm = TRUE) # normalize by median weight
  
  return(val)
}

gradient_t <- function(X, Yhat, wt, lambda, vm, vs, vt,
                       normalize = TRUE) {
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(dim(X)[1] != length(vm))
    stop("Dimensions of X and vm do not agree.")
  if(dim(X)[2] != length(vs))
    stop("Dimensions of X and vs do not agree.")
  if(dim(X)[3] != length(vt))
    stop("Dimensions of X and vt do not agree.")
  
  ms <- outer(vm, vs)
  Xsum_m <- apply(X, c(2, 3), sum)
  expY <- exp(Yhat + lambda * outer(vm, outer(vs, vt)))
  expYsum_m <- apply(expY, c(2, 3), sum)
  
  val <- - vapply(seq(1, length(vt)),
                  function(k) {
                    sum((apply(lambda * ms * X[, , k], 2, sum) - 
                           Xsum_m[, k] *
                           apply(expY[, , k] * lambda * ms, 2, sum) /
                           expYsum_m[, k]) / 
                          wt[, k],
                        na.rm = TRUE)
                  },
                  0.0)
  if(normalize)
    val <- val / 
    median(Xsum_m, na.rm = TRUE) / # normalize by median library size
    sum(!is.na(Xsum_m)) * # normalize by sample size
    median(wt, na.rm = TRUE) # normalize by median weight
  
  return(val)
}

gradient_lambda <- function(X, Yhat, wt, lambda, vm, vs, vt,
                            normalize = TRUE) {
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(dim(X)[1] != length(vm))
    stop("Dimensions of X and vm do not agree.")
  if(dim(X)[2] != length(vs))
    stop("Dimensions of X and vs do not agree.")
  if(dim(X)[3] != length(vt))
    stop("Dimensions of X and vt do not agree.")
  
  mst <- outer(vm, outer(vs, vt))
  Xsum_m <- apply(X, c(2, 3), sum)
  expY <- exp(Yhat + lambda * mst)
  expYsum_m <- apply(expY, c(2, 3), sum)
  
  
  val <- - sum((apply(X * mst, c(2, 3), sum) - 
                  Xsum_m *
                  apply(expY * mst, c(2, 3), sum) /
                  expYsum_m) /
                 wt,
               na.rm = TRUE)
  if(normalize)
    val <- val / 
    median(Xsum_m, na.rm = TRUE) / # normalize by median library size
    sum(!is.na(Xsum_m)) * # normalize by sample size
    median(wt, na.rm = TRUE) # normalize by median weight
  
  return(val)
}
