//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include "miceFast.h"

//weighted linear regression

arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X,arma::colvec &w, arma::mat &X1) {

  int n = X.n_rows, k = X.n_cols;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(n,k);

  for(int h=0;h<k;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2, y2);

  return X1*coef;

}

//weighted linear regression - noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(N,C);

  for(int h=0;h<C;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2, y2);

  arma::colvec res = y - X*coef;

  double sigma = sqrt(arma::as_scalar(arma::trans(res)*res)/(N - C - 1));

  arma::vec noise(N_NA);
  noise.randn();

  return X1 * coef + noise * sigma;

}

//weighted linear regression - bayes

arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(N,C);

  for(int h=0;h<C;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2, y2);

  arma::colvec res = y - X*coef;

  double df = N-C ;

  arma::vec chi2 = Rcpp::rchisq(1, df);

  double sigma_b = sqrt(arma::as_scalar(arma::trans(res)*res)/arma::as_scalar(chi2));

  arma::vec noise(C);
  noise.randn();

  arma::mat XX_inv = arma::inv(arma::trans(X)*X) ;

  arma::colvec beta_b = coef + arma::trans(arma::chol(XX_inv)) * noise * sigma_b;

  arma::vec noise2(N_NA);
  noise2.randn();

  return X1 * coef + noise2 * sigma_b;

}

//linear regression - bayes

arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec coef = arma::solve(X,y);//,arma::solve_opts::fast);//arma::inv(X.t()*X)*X.t()*y;

  arma::colvec res = y - X*coef;

  double df = N-C ;

  arma::vec chi2 = Rcpp::rchisq(1, df);

  double sigma_b = sqrt(arma::as_scalar(arma::trans(res)*res)/arma::as_scalar(chi2));

  arma::vec noise(C);
  noise.randn();

  arma::mat XX_inv = arma::inv(arma::trans(X)*X) ;

  arma::colvec beta_b = coef + arma::trans(arma::chol(XX_inv)) * noise * sigma_b;

  arma::vec noise2(N_NA);
  noise2.randn();

  return X1 * coef + noise2 * sigma_b;

}

//Linear regression with noise
arma::colvec fastLm_noise(arma::colvec &y,arma::mat &X, arma::mat &X1) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec coef = arma::solve(X,y);//,arma::solve_opts::fast);//arma::inv(X.t()*X)*X.t()*y;

  arma::colvec res = y - X*coef;

  double sigma = sqrt(arma::as_scalar(arma::trans(res)*res)/(N - C - 1));

  arma::vec noise(N_NA);
  noise.randn();

  return X1 * coef + noise * sigma;

}

//Simple linear regression

arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1) {

  arma::colvec coef = arma::solve(X, y);//,arma::solve_opts::fast);//arma::inv(X.t()*X)*X.t()*y;

  return X1*coef;

}

//LDA prediction model

arma::colvec fastLda( arma::colvec &y,  arma::mat &X, arma::mat &X1) {

  arma::uvec vars = arma::find(arma::var(X)>0);

  X = X.cols(vars);

  double tol = 1e-6;

  int N = X.n_rows; int C = X.n_cols;

  arma::vec un = arma::unique(y);

  int group = un.n_elem;

  if(group==1 || group >15){
    Rcpp::stop("minimum 2 and maximum 15 categories");
  }

  arma::vec counts = arma::conv_to<arma::vec>::from(arma::hist(y,un));

  arma::vec prior = counts/ (double) N;

  arma::mat group_means(group,C);

  for(int i=0;i<group;i++){

    arma::uvec index = arma::find(y == un(i));
    group_means.row(i) = arma::mean(X.rows(index),0);

  }

  arma::mat group_means_mat(N,C);

  for(int i=0;i<N;i++) group_means_mat.row(i) = group_means.rows(arma::find(un == y(i)));

  arma::mat Sw = arma::sqrt(arma::cov(X - group_means_mat));

  arma::vec scaling0 = 1/Sw.diag();

  arma::mat scaling = arma::diagmat(scaling0);

  double fac =  1/(double)N;

  arma::mat X0 = sqrt(fac) * (X - group_means_mat) * scaling;

  arma::mat input = X0.t()*X0;
  arma::mat U;
  arma::mat V;
  arma::vec s;

  arma::svd( U, s, V, input);

  arma::uvec proper = arma::find(s > tol);

  int rank = proper.n_elem;

  if(rank == 0){
    Rcpp::stop("rank = 0: variables are numerically constant");
  }

  if(rank < C){
    Rcpp::warning("variables are collinear");
  }

  scaling = scaling * U.cols(proper)  * arma::diagmat(1/s.elem(proper));

  arma::vec xb = group_means.t()*prior;

  fac = 1/(double) C;

  X0 = arma::conv_to<arma::rowvec>::from(sqrt(((double) N * prior)*fac))* (group_means.each_row() - xb.t())  * scaling;

  arma::svd( U, s, V, X0.t()*X0 );

  proper = arma::find(s > tol * s(0));

  rank = proper.n_elem;

  if(rank == 0){
    Rcpp::stop("group means are numerically identical");
  }

  scaling = scaling * U.cols(proper);

  arma::mat X_means = arma::mean(X,0);

  arma::vec dist_base = 0.5 * arma::mean(arma::pow(group_means * scaling,2),1) - arma::log(prior);

  N = X1.n_rows;

  //Prediction

  X1 = X1.cols(vars);

  arma::mat dist_base_mat =  arma::repmat(arma::trans(dist_base),N,1);

  arma::mat dist_raw = dist_base_mat - X1 * scaling * arma::trans(group_means * scaling);

  arma::mat dist_scale = arma::exp(-(dist_raw - arma::repmat(arma::min(dist_raw,1),1,group)));

  arma::colvec pred(N);

  for(int i =0;i<N;i++){
    int pos = arma::index_max(dist_scale.row(i));
    pred(i) = un(pos);
  }

  return pred;


}
