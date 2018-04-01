
/*-------------------------------------------------------------------------------

 This file is part of miceFast.

 miceFast is free software: you can redistribute it and/or modify

 it under the terms of the GNU General Public License as published by

 the Free Software Foundation, either version 3 of the License, or

 (at your option) any later version.


 miceFast is distributed in the hope that it will be useful,

 but WITHOUT ANY WARRANTY; without even the implied warranty of

 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the

 GNU General Public License for more details.


 You should have received a copy of the GNU General Public License

 along with miceFast. If not, see <http://www.gnu.org/licenses/>.


 Written by:

 Maciej Nasinski

#-------------------------------------------------------------------------------*/
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include "miceFast.h"

arma::uvec complete_cases_mat(arma::mat &x) {

    unsigned int N = x.n_rows;

    arma::uvec x_out(N,arma::fill::ones);

    for(unsigned int a=0;a<N;a++){

      bool na = x.row(a).has_nan();

      if(na){
        x_out(a) = 0;
      }

    }

    return x_out;

  }

arma::uvec complete_cases_vec(arma::colvec &y) {

  unsigned int N = y.n_elem;

  arma::uvec out(N,arma::fill::ones);

  for(unsigned int i = 0; i<N; i++){

    if(y.row(i).has_nan()){out(i)=0;}

  }

  return out;
}

arma::uvec histFast(arma::uvec &gg){
  unsigned int n = gg.n_elem;
  arma::uvec un = arma::unique(gg);
  unsigned int group = un.n_elem;
  arma::uvec h(group,arma::fill::zeros);
  unsigned int gg_prev,gg_curr;
  gg_prev = gg(0);
  unsigned int iter = 0 ;
  h(iter) = 0;
  for(unsigned int i=0;i<n;i++){
    gg_curr = gg(i);
    if(gg_prev!=gg_curr){
      iter++;
    }
    h(iter) = h(iter) + 1;
    gg_prev = gg_curr;
  }
  return h;
}


bool different_y_and_x(int posit_y, arma::uvec posit_x){

  arma::uvec different_ind = arma::find(posit_x==posit_y);

  bool different_vars = different_ind.n_elem == 0;

  return different_vars;
}

bool different_x(arma::uvec posit_x){

  arma::uvec uni_x = arma::unique(posit_x);

  bool different_vars = uni_x.n_elem == posit_x.n_elem;

  return different_vars;
}
