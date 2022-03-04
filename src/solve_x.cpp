#include <RcppArmadillo.h>
#include "logLik.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List solve_x_vm(const double& L_old,
                      const double& gamma,
                      const arma::vec& gradient,
                      const arma::cube& X,
                      const arma::cube& Yhat,
                      const arma::mat& wt,
                      const double& lambda,
                      const arma::vec& vm,
                      const arma::vec& vs,
                      const arma::vec& vt,
                      const unsigned int& max_iter = 1000) {
  if(gradient.n_elem != vm.n_elem) {
    throw std::invalid_argument("Dimensions of x and gradient don't agree!");
  }
  
  arma::vec x_old = vm;
  double f_old = negLogLik_cpp(X, Yhat, wt, lambda, x_old, vs, vt);
  double L_new = L_old, f_new, diff;
  arma::vec i_step = gradient / L_new;
  double i_term = arma::accu(arma::pow(gradient, 2)) / L_new / 2;
  arma::vec x_new;
  unsigned int i_iter;
  
  for(i_iter = 0; i_iter <= max_iter; i_iter++) {
    x_new = x_old - i_step;
    
    f_new = negLogLik_cpp(X, Yhat, wt, lambda, x_new, vs, vt);
    diff = f_new - f_old + i_term;
    
    if(diff <= 0) {
      break;
    }
    
    L_new *= gamma;
    i_step /= gamma;
    i_term /= gamma;
  }
  
  return Rcpp::List::create(Rcpp::Named("x_new") = x_new,
                            Rcpp::Named("L_new") = L_new,
                            Rcpp::Named("f_new") = f_new,
                            Rcpp::Named("n_iter") = i_iter);
}

// [[Rcpp::export]]
Rcpp::List test_cpp(const arma::cube& X) {
  
  unsigned int n_m = X.n_rows, n_s = X.n_cols, 
    n_t = X.n_slices;
  arma::mat Xsum_m(n_s, n_t);
  Xsum_m = arma::sum(X, 0);
  
  return Rcpp::List::create(Rcpp::Named("test1")=Xsum_m);
}
  
// solve_x <- function(L_old, gamma, gradient, f, x_old,
//                     max_iter = 1000) {
//   if(length(x_old) != length(gradient))
//     stop("Dimensions of x_old and gradient don't agree!")
//     
//     f_old <- f(x_old)
//     
//     i_try <- 0
//   while(TRUE) {
//     L_new <- L_old * gamma^i_try
//     x_new <- x_old - gradient / L_new
//     
//     f_new <- f(x_new)
//     diff <- f_new - f_old + sum(gradient^2) / L_new / 2
//     
//     if(diff <= 0 | i_try > max_iter)
//       break
//       
//       i_try <- i_try + 1
//   }
//   
//   return(list(x_new = x_new,
//               L_new = L_new,
//               f_new = f_new))
// }