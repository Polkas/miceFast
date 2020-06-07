//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include "miceFast.h"
#include <stdlib.h>     /* rand */
#include <algorithm>
using namespace std;
using namespace Rcpp;

#include <algorithm>
using namespace std;
using namespace Rcpp;
using namespace arma;


#define UNUSED(expr) (void)(expr)

//weighted linear regression
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X,arma::colvec &w, arma::mat &X1,int times) {

  UNUSED(times);

  int n = X.n_rows, k = X.n_cols;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(n,k);

  for(int h=0;h<k;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2.t()*X2, X2.t()*y2);

  return X1*coef;

}

//weighted linear regression - noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1,int times) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(N,C);

  for(int h=0;h<C;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2.t()*X2, X2.t()*y2);

  arma::colvec res = y - X*coef;

  double sigma = sqrt(arma::as_scalar(arma::trans(res)*res)/(N - C - 1));

  arma::colvec pred_sum(N_NA,arma::fill::zeros);

  for(int i=0;i<times;i++){

    arma::vec noise(N_NA);
    noise.randn();

    pred_sum = pred_sum + (X1 * coef + noise * sigma);
  }

  return pred_sum/(double) times;

}

//weighted linear regression - bayes

arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1,int times) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(N,C);

  for(int h=0;h<C;h++){
    X2.col(h)=wq%X.col(h);
  }

  arma::colvec coef = arma::solve(X2.t()*X2, X2.t()*y2);

  arma::colvec res = y - X*coef;

  double df = N-C ;

  //arma::mat XX_inv = arma::inv(arma::trans(X)*X) ;

  double res2 = arma::as_scalar(arma::trans(res)*res);

  arma::colvec pred_sum(N_NA,arma::fill::zeros);

  for(int i=0;i<times;i++){

    double chi2 = Rcpp::as<double>(Rcpp::rchisq(1, df));

    double sigma_b = sqrt(res2/chi2);

    arma::vec noise2(N_NA);
    noise2.randn();

    pred_sum = pred_sum +( X1 * coef + noise2 * sigma_b);
  }

  return pred_sum/(double) times;

}

//linear regression - bayes

arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1,int times) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec coef =  arma::solve(X.t()*X, X.t()*y);

  arma::colvec res = y - X*coef;

  double df = N-C;

  //arma::mat XX_inv = arma::inv(arma::trans(X)*X) ;

  double res2 = arma::as_scalar(arma::trans(res)*res);

  arma::colvec pred_sum(N_NA,arma::fill::zeros);

  for(int i=0;i<times;i++){

    double chi2 = Rcpp::as<double>(Rcpp::rchisq(1, df));

    double sigma_b = sqrt(res2/chi2);

    arma::vec noise2(N_NA);
    noise2.randn();

    pred_sum = pred_sum + (X1 * coef + noise2 * sigma_b);
  }

  return pred_sum/(double) times;

}

//Linear regression with noise
arma::colvec fastLm_noise(arma::colvec &y,arma::mat &X, arma::mat &X1,int times) {

  int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;

  arma::colvec coef =  arma::solve(X.t()*X, X.t()*y);

  arma::colvec res = y - X*coef;

  double sigma = sqrt(arma::as_scalar(arma::trans(res)*res)/(N - C - 1));

  arma::colvec pred_sum(N_NA,arma::fill::zeros);

  for(int i=0;i<times;i++){

    arma::vec noise(N_NA);
    noise.randn();

    pred_sum = pred_sum + (X1 * coef + noise * sigma);
  }

  return pred_sum/(double) times;

}

//Simple linear regression

arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1,int times) {

  UNUSED(times);

  arma::colvec coef = arma::solve(X.t()*X, X.t()*y);

  return X1*coef;

}

//LDA prediction model

arma::colvec fastLda( arma::colvec &y,  arma::mat &X, arma::mat &X1, int times) {

  UNUSED(times);

  arma::uvec vars = arma::find(arma::var(X)>0);

  arma::mat X_vol = X.cols(vars);

  arma::mat X_means = arma::mean(X_vol,0);

  double tol = 1e-6;

  int N = X_vol.n_rows; int C = X_vol.n_cols;

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
    group_means.row(i) = arma::mean(X_vol.rows(index),0);

  }

  arma::mat group_means_mat(N,C);

  for(int i=0;i<N;i++) group_means_mat.row(i) = group_means.rows(arma::find(un == y(i)));

  arma::mat Sw = arma::sqrt(arma::cov(X_vol - group_means_mat));

  arma::vec scaling0 = 1/Sw.diag();

  arma::mat scaling = arma::diagmat(scaling0);

  double fac =  1/(double)N;

  arma::mat X0 = sqrt(fac) * (X_vol - group_means_mat) * scaling;

  X_vol.clear();

  arma::mat input = X0.t()*X0;
  X0.clear();
  arma::mat U;
  arma::mat V;
  arma::vec s2;

  arma::svd( U, s2, V, input);

  arma::vec s = arma::sqrt(s2);

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

  arma::mat X_s = arma::conv_to<arma::rowvec>::from(sqrt(((double) N * prior)*fac))* (group_means.each_row() - xb.t())  * scaling;

  arma::svd( U, s2, V, X_s.t()*X_s );
  X_s.clear();

  s = arma::sqrt(s2);

  proper = arma::find(s > tol * s(0));

  rank = proper.n_elem;

  if(rank == 0){
    Rcpp::stop("group means are numerically identical");
  }

  scaling = scaling * U.cols(proper);

  arma::vec dist_base = 0.5 * arma::mean(arma::pow(group_means * scaling,2),1) - arma::log(prior);

  N = X1.n_rows;

  //Prediction

  arma::mat X1_vol = X1.cols(vars);

  arma::mat dist_base_mat =  arma::repmat(arma::trans(dist_base),N,1);

  arma::mat dist_raw = dist_base_mat - X1_vol * scaling * arma::trans(group_means * scaling);

  dist_raw.each_col() -= arma::min(dist_raw,1);

  dist_raw.transform( [](double val) { return (-std::exp(val)); } );

  arma::colvec pred(N);

  for(int i =0;i<N;i++){
    int pos = arma::index_max(dist_raw.row(i));
    pred(i) = un(pos);
  }

  return pred;


}


//QDA prediction model

//PCA prediction model

//Ridge prediction model

//PMM

int findCrossOver(arma::colvec arr, double low, double high, double x)

{
  if (arr[high] <= x)
    return high;
  if (arr[low] > x)
    return low;


  int mid = (low + high)/2;

  if (arr[mid] <= x && arr[mid+1] > x)
    return mid;

  if(arr[mid] < x)
    return findCrossOver(arr, mid+1, high, x);

  return findCrossOver(arr, low, mid - 1, x);
}


arma::colvec Kclosestrand(arma::colvec arr, double x, int k)
{
  int n = arr.n_rows;


  int l = findCrossOver(arr, 0, n-1, x);
  int r = l;
  int count = 0;
  arma::colvec resus(k);


  if (arr[l] == x) l--;

  while (l >= 0 && r < n && count < k)
  {
    if (x - arr[l] < arr[r] - x)
      resus[count] = arr[l--];
    else
      resus[count] = arr[r++];
    count++;
  }

  while (count < k && l >= 0)
    resus[count] = arr[l--], count++;

  while (count < k && r < n)
    resus[count] = arr[r++], count++;

  return resus;

}



//' Finding in random manner one of the k closets points in a certain vector for each value in a second vector
//'
//' @description this function using pre-sorting of a \code{y} and the by the binary search the one of the k closest value for each miss is returned.
//'
//' @param y numeric vector values to be look up
//' @param miss numeric vector a values to be look for
//' @param k integer a number of values which should be taken into account during sampling one of the k closest point
//'
//' @return a numeric vector

//'
//' @name neibo
//'
//' @export
// [[Rcpp::export]]
arma::colvec neibo(arma::colvec y, arma::colvec miss, int k) {

  int n_y = y.n_rows;

  k = (k <= n_y) ? k : n_y;
  k = (k >= 1) ? k : 1;

  arma::colvec y_new = arma::sort(y);

  int n_miss = miss.n_rows;

  arma::colvec resus(n_miss);
  arma::colvec resus_temp(k);

  arma::vec which = floor(runif(n_miss, 0, k));


    int subs;

  for(int i=0; i<n_miss ;i++){
    double mm = miss[i];
    resus_temp = Kclosestrand(y_new,mm,k);
    subs = (int) which[i];
    resus[i] = resus_temp[subs];
  }

  return resus ;

}



arma::uvec Kclosestrand_index(arma::colvec arr, double x, int k)
{
  int n = arr.size();

  int l = findCrossOver(arr, 0, n-1, x);
  int r = l;
  int count = 0;
  arma::uvec resus(k);

  if (arr[l] == x) l--;

  while (l >= 0 && r < n && count < k)
  {
    if (x - arr[l] < arr[r] - x)
      resus[count] = l--;
    else
      resus[count] = r++;
    count++;
  }

  while (count < k && l >= 0)
    resus[count] = l--, count++;

  while (count < k && r < n)
    resus[count] = r++, count++;

  return resus;

}



//' Finding in random manner one of the k closets points in a certain vector for each value in a second vector
//'
//' @description this function using pre-sorting of a \code{y} and the by the binary search the one of the k closest value for each miss is returned.
//'
//' @param y numeric vector values to be look up
//' @param miss numeric vector a values to be look for
//' @param k integer a number of values which should be taken into account during sampling one of the k closest point
//'
//' @return a numeric vector
//'
//' @name neibo_index
//'
//' @export
// [[Rcpp::export]]
arma::uvec neibo_index(arma::colvec y, arma::colvec miss, int k) {
  int n_y = y.n_rows;

  k = (k <= n_y) ? k : n_y;
  k = (k >= 1) ? k : 1;

  arma::colvec y_new = arma::sort(y);

  int n_miss = miss.n_rows;

  arma::uvec resus(n_miss);
  arma::uvec resus_temp(k);

  arma::vec which = floor(runif(n_miss, 0, k));


  int subs;

  for(int i=0; i<n_miss ;i++){
    double mm = miss[i];
    resus_temp = Kclosestrand_index(y_new,mm,k);
    subs = (int) which[i];
    resus[i] = resus_temp[subs];
  }

  return resus ;

}


// [[Rcpp::export]]
arma::colvec pmm_weighted_neibo( arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1,int k) {

  int N = X.n_rows; int C = X.n_cols;

  arma::colvec wq=sqrt(w);
  arma::colvec y2 = wq%y;
  arma::mat X2(N,C);

  for(int h=0;h<C;h++){
    X2.col(h)=wq%X.col(h);
  }


    arma::mat xtx = arma::mat( arma::trans(X2) * X2 );
    for (int ii=0;ii<k;ii++){
        xtx(ii,ii)=xtx(ii,ii)+ridge;
    }
    arma::mat xinv = arma::inv( xtx );
    arma::mat coef2 = arma::mat( xinv * arma::trans(X2) * y );
    arma::colvec resid = arma::mat( y - X*coef2 );
    double sig2 = arma::as_scalar( arma::trans(resid)*resid/(N-C) );

    arma::mat vcoef = arma::mat( sig2 * xinv );

    arma::colvec coefu(k,fill::randn);

    arma::mat coef = arma::mat( coef2 + arma::chol(vcoef) * coefu );


    arma::colvec ypred_mis = X1 * coef;
    arma::colvec y_full = arma::sort(y);

    arma::colvec yimp = neibo(y_full,ypred_mis,k);

    return yimp;
}

arma::colvec pmm_neibo( arma::colvec &y, arma::mat &X,arma::mat &X1,int k) {

  int N = X.n_rows; int C = X.n_cols;


    arma::mat xtx = arma::mat( arma::trans(X) * X );
    for (int ii=0;ii<k;ii++){
        xtx(ii,ii)=xtx(ii,ii)+ridge;
    }
    arma::mat xinv = arma::inv( xtx );
    arma::mat coef2 = arma::mat( xinv * arma::trans(X) * y );
    arma::colvec resid = arma::mat( y - X*coef2 );
    double sig2 = arma::as_scalar( arma::trans(resid)*resid/(N-C) );

        arma::mat vcoef = arma::mat( sig2 * xinv );

    arma::colvec coefu(k,fill::randn);

    arma::mat coef = arma::mat( coef2 + arma::chol(vcoef) * coefu );

   arma::colvec ypred_mis = X1 * coef;
    arma::colvec y_full = arma::sort(y);


    arma::colvec yimp = neibo(y_full,ypred_mis,k);

    return yimp;
}

static R_CallMethodDef callMethods[]  = {
  {"neibo", (DL_FUNC) &neibo, 3},
  {"neibo_index", (DL_FUNC) &neibo_index, 3},
  {"pmm_neibo", (DL_FUNC) &pmm_neibo, 4},
  {"pmm_weighted_neibo", (DL_FUNC) &pmm_weighted_neibo, 5},
  {NULL, NULL, 0}
};
