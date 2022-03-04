#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector gradient_m_cpp(const arma::cube& X,
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
  
  arma::mat st = vs * vt.t();
  arma::cube mst(arma::size(X));
  for(unsigned int i_m = 0; i_m < n_m; i_m++) {
    mst.row(i_m) = vm(i_m) * st;
  }
  arma::mat Xsum_m = arma::sum(X, 0);
  arma::cube Yhat_new = Yhat + lambda * mst;
  arma::cube expY = arma::exp(Yhat_new);
  arma::mat expYsum_m = arma::sum(expY, 0);
  arma::mat term_prod = Xsum_m / expYsum_m;
  
  arma::vec val(vm.n_elem);
  for(unsigned int i_m = 0; i_m < vm.n_elem; i_m++) {
    arma::mat i_m_X = X.row(i_m);
    arma::mat i_m_expY = expY.row(i_m);
    arma::mat i_mat_val = 
      lambda * st % (i_m_X - term_prod % i_m_expY) /
      wt;
    val(i_m) = -arma::accu(i_mat_val.replace(arma::datum::nan, 0));
  }
  
  if(normalize) {
    arma::uvec non_missing = arma::find_finite(Xsum_m);
    val = val / arma::median(Xsum_m.elem(non_missing)) /
      non_missing.n_elem *
        arma::median(wt.elem(non_missing));  
  }
  
  return Rcpp::NumericVector(val.begin(),val.end());
}


// [[Rcpp::export]]
Rcpp::NumericVector gradient_s_cpp(const arma::cube& X,
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
  
  arma::mat mt = vm * vt.t();
  arma::cube mst(arma::size(X));
  for(unsigned int i_s = 0; i_s < n_s; i_s++) {
    mst.col(i_s) = vs(i_s) * mt;
  }
  arma::mat Xsum_m = arma::sum(X, 0);
  arma::cube Yhat_new = Yhat + lambda * mst;
  arma::cube expY = arma::exp(Yhat_new);
  arma::mat expYsum_m = arma::sum(expY, 0);
  
  arma::vec val(vs.n_elem);
  for(unsigned int i_s = 0; i_s < vs.n_elem; i_s++) {
    arma::mat i_s_X = X.col(i_s);
    arma::mat i_s_expY = expY.col(i_s);
    arma::rowvec i_vec_val = (
      arma::sum(lambda * mt % i_s_X, 0) -
        Xsum_m.row(i_s) %
        arma::sum(lambda * i_s_expY % mt, 0) / expYsum_m.row(i_s)
    ) / wt.row(i_s);
    val(i_s) = -arma::accu(i_vec_val.replace(arma::datum::nan, 0));
  }
  
  if(normalize) {
    arma::uvec non_missing = arma::find_finite(Xsum_m);
    val = val / arma::median(Xsum_m.elem(non_missing)) /
      non_missing.n_elem *
        arma::median(wt.elem(non_missing));
  }

  return Rcpp::NumericVector(val.begin(),val.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector gradient_t_cpp(const arma::cube& X,
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
  
  arma::mat ms = vm * vs.t();
  arma::cube mst(arma::size(X));
  for(unsigned int i_t = 0; i_t < n_t; i_t++) {
    mst.slice(i_t) = vt(i_t) * ms;
  }
  arma::mat Xsum_m = arma::sum(X, 0);
  arma::cube Yhat_new = Yhat + lambda * mst;
  arma::cube expY = arma::exp(Yhat_new);
  arma::mat expYsum_m = arma::sum(expY, 0);
  
  arma::vec val(vt.n_elem);
  for(unsigned int i_t = 0; i_t < vt.n_elem; i_t++) {
    arma::mat i_t_X = X.slice(i_t);
    arma::mat i_t_expY = expY.slice(i_t);
    arma::rowvec i_vec_val = (
      arma::sum(lambda * ms % i_t_X, 0) -
        Xsum_m.col(i_t).t() %
        arma::sum(lambda * i_t_expY % ms, 0) / expYsum_m.col(i_t).t()
    ) / wt.col(i_t).t();
    val(i_t) = -arma::accu(i_vec_val.replace(arma::datum::nan, 0));
  }
  
  if(normalize) {
    arma::uvec non_missing = arma::find_finite(Xsum_m);
    val = val / arma::median(Xsum_m.elem(non_missing)) /
      non_missing.n_elem *
        arma::median(wt.elem(non_missing));
  }
  
  return Rcpp::NumericVector(val.begin(),val.end());
}

// [[Rcpp::export]]
double gradient_lambda_cpp(const arma::cube& X,
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
  
  arma::mat X_mst_sum = arma::sum(X % mst, 0);
  arma::mat expY_mst_sum = arma::sum(expY % mst, 0);
  arma::mat val_sum = (X_mst_sum -
    Xsum_m % expY_mst_sum / expYsum_m) / 
    wt;
  double val = -arma::accu(val_sum.replace(arma::datum::nan, 0));
  
  if(normalize) {
    arma::uvec non_missing = arma::find_finite(Xsum_m);
    val = val / arma::median(Xsum_m.elem(non_missing)) /
      non_missing.n_elem *
        arma::median(wt.elem(non_missing));
  }
  
  return val;
}