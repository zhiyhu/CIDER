#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate cosine simialrity
//' @description Refer to https://gist.github.com/bobthecat/2903031
//' @param Xr a matrix
//' @export
// [[Rcpp::export]]
NumericMatrix cosine_similarity(NumericMatrix Xr) {
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);
  arma::mat Y = arma::trans(X) * X; // matrix product
  arma::mat res = (1 - Y / (arma::sqrt(arma::diagvec(Y)) * arma::trans(arma::sqrt(arma::diagvec(Y)))));
  return Rcpp::wrap(res);
}