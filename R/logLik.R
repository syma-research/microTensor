# negative (quasi-)log likelihood for the microTensor model
negLogLik <- function(X, Yhat, wt, lambda, vm, vs, vt, normalize = TRUE) {
  if(!all(dim(X) == dim(Yhat)))
    stop("Dimensions of X and Yhat do not agree.")
  if(!all(dim(X)[c(2, 3)] == dim(wt)))
    stop("Dimensions of X and wt do not agree.")
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

  val <- sum((-apply((Yhat + lambda * mst) * X, c(2 ,3), sum) + 
                Xsum_m * log(expYsum_m)) / wt,
             na.rm = TRUE)
  if(normalize)
    val <- val / 
    median(Xsum_m, na.rm = TRUE) / # normalize by median library size
    sum(!is.na(Xsum_m)) * # normalize by sample size
    median(wt, na.rm = TRUE) # normalize by median weight
  
  return(val)
}
