//
//  getP.cpp
//
//
//  Created by Thijs Janzen on 11/04/2019.
//

#include <Rcpp.h>
#include <cmath>

/*
// from phangorn:
void getP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double res; //tmp[m],
    double *tmp;
    tmp = (double *) malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++) res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
    free(tmp);  // ausserhalb
}
// 64 * 3 + 4     16 + 64 * 2
*/

// [[Rcpp::export]]
Rcpp::NumericMatrix getPM_rcpp(Rcpp::List eig,
                               int el, //double branch_length,
                               double w) //double rate) {
{
  Rcpp::NumericVector eva   = Rcpp::as<Rcpp::NumericVector>(eig["values"]);
  Rcpp::NumericVector ev    = Rcpp::as<Rcpp::NumericMatrix>(eig["vectors"]);
  Rcpp::NumericVector evei  = Rcpp::as<Rcpp::NumericMatrix>(eig["inv"]);
  Rcpp::NumericMatrix P(4, 4);
  if(el == 0 || w == 0) {
      return P;
  }

  Rcpp::NumericVector tmp(4);
  for(int i = 0; i < 4; i++) tmp[i] = std::exp(1.0 * eva[i] * w * el);

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


