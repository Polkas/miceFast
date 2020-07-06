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

// arma::uvec histFast(arma::uvec &gg){
//   unsigned int n = gg.n_elem;
//   arma::uvec un = arma::unique(gg);
//   unsigned int group = un.n_elem;
//   arma::uvec h(group,arma::fill::zeros);
//   unsigned int gg_prev,gg_curr;
//   gg_prev = gg(0);
//   unsigned int iter = 0 ;
//   h(iter) = 0;
//   for(unsigned int i=0;i<n;i++){
//     gg_curr = gg(i);
//     if(gg_prev!=gg_curr){
//       iter++;
//     }
//     h(iter) = h(iter) + 1;
//     gg_prev = gg_curr;
//   }
//   return h;
// }

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

arma::mat sym(arma::mat x){

  return (x + arma::trans(x))/2;

}
