#ifndef MODSTRING_H
#define MODSTRING_H

#include <RcppArmadillo.h>
double negLogLik_cpp(const arma::cube& X,
                     const arma::cube& Yhat,
                     const arma::mat& wt,
                     const double& lambda,
                     const arma::vec& vm,
                     const arma::vec& vs,
                     const arma::vec& vt,
                     const bool& normalize = true);

#endif