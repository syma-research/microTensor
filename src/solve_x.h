#ifndef MODSTRING_H
#define MODSTRING_H

#include <RcppArmadillo.h>
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
                      const unsigned int& max_iter = 1000);

#endif