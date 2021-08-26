# coordinate gradient descent for each parameter

update_m <- function(vm_old, X, Yhat, wt, m_proj,
                     lambda, vs, vt, 
                     gamma = 1.25, L_old = 1) {
  f <- function(x) 
    negLogLik(X = X, Yhat = Yhat, wt = wt, 
              lambda = lambda, vm = x, vs = vs, vt = vt,
              normalize = TRUE)
  
  # update
  gradient <- gradient_m(X = X, Yhat = Yhat, wt = wt, 
                         lambda = lambda, vm = vm_old, vs = vs, vt = vt,
                         normalize = TRUE)
  l_update <- solve_x(L_old = L_old,
                      gamma = gamma,
                      gradient = gradient, 
                      f = f,
                      x_old = vm_old)
  vm_new <- l_update$x_new
  
  # projection
  vm_new <- as.vector(m_proj %*% vm_new)
  
  return(list(vm = vm_new,
              L = l_update$L_new,
              f = f(vm_new)))
}

update_s <- function(vs_old, X, Yhat, wt, 
                     lambda, vm, vt, 
                     gamma = 1.25, L_old = 1) {
  f <- function(x) 
    negLogLik(X = X, Yhat = Yhat, wt = wt, 
              lambda = lambda, vm = vm, vs = x, vt = vt,
              normalize = TRUE)
  
  # update
  gradient <- gradient_s(X = X, Yhat = Yhat, wt = wt, 
                         lambda = lambda, vm = vm, vs = vs_old, vt = vt,
                         normalize = TRUE)
  l_update <- solve_x(L_old = L_old,
                      gamma = gamma,
                      gradient = gradient, 
                      f = f,
                      x_old = vs_old)
  vs_new <- l_update$x_new
  
  return(list(vs = vs_new,
              L = l_update$L_new,
              f = f(vs_new)))
}

update_t <- function(vt_old, X, Yhat, wt, nn_t = FALSE,
                     lambda, vm, vs, 
                     gamma = 1.25, L_old = 1) {
  f <- function(x) 
    negLogLik(X = X, Yhat = Yhat, wt = wt, 
              lambda = lambda, vm = vm, vs = vs, vt = x,
              normalize = TRUE)
  
  # update
  gradient <- gradient_t(X = X, Yhat = Yhat, wt = wt, 
                         lambda = lambda, vm = vm, vs = vs, vt = vt_old,
                         normalize = TRUE)
  l_update <- solve_x(L_old = L_old,
                      gamma = gamma,
                      gradient = gradient,
                      f = f,
                      x_old = vt_old)
  vt_new <- l_update$x_new
  
  # projection
  if(nn_t)
    vt_new[vt_new < 0] <- 0
  
  return(list(vt = vt_new,
              L = l_update$L_new,
              f = f(vt_new)))
}

update_lambda <- function(lambda_old, X, Yhat, wt, vm, vs, vt,
                          gamma = 1.25, L_old = 1) {
  f <- function(x) 
    negLogLik(X = X, Yhat = Yhat, wt = wt, 
              lambda = x, vm = vm, vs = vs, vt = vt,
              normalize = TRUE)
  
  # update
  gradient <- gradient_lambda(X = X, Yhat = Yhat, wt = wt, 
                              lambda = lambda_old, vm = vm, vs = vs, vt = vt,
                              normalize = TRUE)
  l_update <- solve_x(L_old = L_old,
                      gamma = gamma,
                      gradient = gradient, 
                      f = f,
                      x_old = lambda_old)
  lambda_new <- l_update$x_new
  
  return(list(lambda = lambda_new,
              L = l_update$L_new,
              f = f(lambda_new)))
}

solve_x <- function(L_old, gamma, gradient, f, x_old,
                    max_iter = 1000) {
  if(length(x_old) != length(gradient))
    stop("Dimensions of x_old and gradient don't agree!")
  
  f_old <- f(x_old)
  
  i_try <- 0
  while(TRUE) {
    L_new <- L_old * gamma^i_try
    x_new <- x_old - gradient / L_new
    
    f_new <- f(x_new)
    diff <- f_new - f_old + sum(gradient^2) / L_new / 2
    
    if(diff <= 0 | i_try > max_iter)
        break
    
    i_try <- i_try + 1
  }
  
  return(list(x_new = x_new,
              L_new = L_new,
              f_new = f_new))
}
