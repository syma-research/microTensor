#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double negLogLik_cpp(const arma::cube& X,
                     const arma::cube& Yhat,
                     const arma::mat& wt,
                     const double& lambda,
                     const arma::vec& vm,
                     const arma::vec& vs,
                     const arma::vec& vt,
                     const bool& normalize = true) {
  unsigned int n_m = X.n_rows, n_s = X.n_cols, 
    n_t = X.n_slices;
  if(n_m <= 1 || n_t <= 1 || n_s <= 1) {
    throw std::invalid_argument("Data does not form a tensor!");
  }
  
  if(n_m != Yhat.n_rows ||
     n_s != Yhat.n_cols ||
     n_t != Yhat.n_slices) {
    throw std::invalid_argument("Dimensions of X and Yhat do not agree!");
  }
  if(n_s != wt.n_rows ||
     n_t != wt.n_cols) {
    throw std::invalid_argument("Dimensions of X and wt do not agree!");
  }
  if(n_m != vm.n_rows) {
    throw std::invalid_argument("Dimensions of X and vm do not agree!");
  }
  if(n_s != vs.n_rows) {
    throw std::invalid_argument("Dimensions of X and vs do not agree!");
  }
  if(n_t != vt.n_rows) {
    throw std::invalid_argument("Dimensions of X and vt do not agree!");
  }
  
  
  arma::cube mst(arma::size(X));
  for(unsigned int i_m = 0; i_m < n_m; i_m++) {
    mst.row(i_m) = vm(i_m) * vs * vt.t();
  }
  arma::mat Xsum_m = arma::sum(X, 0);
  arma::cube Yhat_new = Yhat + lambda * mst;
  arma::cube expY = arma::exp(Yhat_new);
  arma::mat expYsum_m = arma::sum(expY, 0);
  
  arma::mat val_sum = -arma::sum(Yhat_new % X, 0);
  val_sum = (val_sum + Xsum_m % arma::log(expYsum_m)) / wt;
  double val = arma::accu(val_sum.replace(arma::datum::nan, 0));
  
  
  if(normalize) {
    arma::uvec non_missing = arma::find_finite(Xsum_m);
    val = val / arma::median(Xsum_m.elem(non_missing)) /
      non_missing.n_elem *
        arma::median(wt.elem(non_missing));  
  }
  
  return val;
}