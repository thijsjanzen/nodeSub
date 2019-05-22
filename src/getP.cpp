//
//  getP.cpp
//
//
//  Created by Thijs Janzen on 11/04/2019.
//

#include <Rcpp.h>
#include <cmath>
// [[Rcpp::export]]
Rcpp::NumericMatrix get_p_m_rcpp(Rcpp::List eig,
                               double branch_length,
                               double rate)
{
  Rcpp::NumericVector eva   = Rcpp::as<Rcpp::NumericVector>(eig["values"]);
  Rcpp::NumericVector ev    = Rcpp::as<Rcpp::NumericMatrix>(eig["vectors"]);
  Rcpp::NumericVector evei  = Rcpp::as<Rcpp::NumericMatrix>(eig["inv"]);

  Rcpp::NumericMatrix local = Rcpp::as<Rcpp::NumericMatrix>(eig["inv"]);
  int dim_size = local.ncol();

  Rcpp::NumericMatrix P(dim_size, dim_size);
  if(branch_length == 0 || rate <= 0) {
      for(int i = 0; i < dim_size; ++i) {
        for(int j = 0; j < dim_size; ++j) {
          if(i != j) P(i,j) = 0;
          if(i == j) P(i,j) = 1;
        }
      }
      return(P);
  }

  Rcpp::NumericVector tmp(dim_size);
  for(int i = 0; i < dim_size; i++) {
    tmp[i] = std::exp(1.0 * eva[i] * rate * branch_length);
  }

  for(int i = 0; i < dim_size; i++){
    for(int j = 0; j < dim_size; j++){
      double res = 0.0;
      for(int h = 0; h < dim_size; h++) {
        res += ev(i, h) * tmp[h] * evei(h, j);
      }
      P(i, j) = res;
    }
  }

  return P;
}
