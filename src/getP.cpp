//
//  getP.cpp
//
//
//  Created by Thijs Janzen on 11/04/2019.
//

#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
Rcpp::NumericMatrix getPM_rcpp(Rcpp::List eig,
                               double branch_length,
                               double rate) {
{
  Rcpp::NumericVector eva   = Rcpp::as<Rcpp::NumericVector>(eig["values"]);
  Rcpp::NumericVector ev    = Rcpp::as<Rcpp::NumericMatrix>(eig["vectors"]);
  Rcpp::NumericVector evei  = Rcpp::as<Rcpp::NumericMatrix>(eig["inv"]);
  Rcpp::NumericMatrix P(4, 4);
  if(branch_length == 0 || rate == 0) {
      return P;
  }

  Rcpp::NumericVector tmp(4);
  for(int i = 0; i < 4; i++) {
    tmp[i] = std::exp(1.0 * eva[i] * rate * branch_length);
  }

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      double res = 0.0;
      for(int h = 0; h < 4; h++) {
        res += ev(i, h) * tmp[h] * evei(h, j);
      }
      P(i, j) = res;
    }
  }

  return P;
}
