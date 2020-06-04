#include <RcppArmadillo.h>
#include "miceFast.h"
#include <math.h>

//functions for getting indexes

arma::uvec get_index_full_R(arma::mat &x,int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index_Full = arma::find((full_y + full_x) == 2);

  return index_Full;
}

arma::uvec get_index_NA_R(arma::mat &x,int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index(x.n_rows,arma::fill::zeros);

  for(unsigned int i=0;i<x.n_rows;i++){

    if( (full_y(i)==0) && (full_x(i)==1) ){index(i) = 1;}

  }

  arma::uvec index_NA = arma::find(index==1);

  return index_NA;
}

arma::colvec impute_raw_R(arma::mat &x,std::string s, int posit_y,arma::uvec posit_x,int times){

  typedef arma::colvec (*pfunc)(arma::colvec&,arma::mat&,arma::mat&,int);
  std::map<std::string, pfunc> funMap = {
    {"lda",fastLda},
    {"lm_pred",fastLm_pred},
    {"lm_noise",fastLm_noise},
    {"lm_bayes",fastLm_bayes}};

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  arma::uvec index_full = get_index_full_R(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA_R(x,posit_y, posit_x);

  arma::colvec pred = x(index_NA,posit_y_uvec);

  if((!(index_NA.n_elem==0)) && ((index_full.n_elem>15 && s=="lda")|| (index_full.n_elem>posit_x.n_elem && s!="lda"))){

  arma::mat X_full = x(index_full,posit_x);
  arma::mat X_NA = x(index_NA,posit_x);
  arma::colvec Y_full = x(index_full,posit_y_uvec);

  pfunc f = funMap[s];
  pred = (*f)(Y_full,X_full,X_NA,times);

  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;

  return Y;
}


arma::colvec imputeW_R(arma::mat &x,std::string s,int posit_y,arma::uvec posit_x,arma::colvec w,int times){

  typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&,int);
  std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
  {"lm_noise",fastLm_weighted_noise},
  {"lm_bayes",fastLm_weighted_bayes}};

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  arma::uvec index_full = get_index_full_R(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA_R(x,posit_y, posit_x);

  if(w.has_nan()){Rcpp::stop("There is NA values for weights");}
  if(arma::any(w<0)){Rcpp::stop("There are negative values for the weights variable");}

  arma::colvec pred = x(index_NA,posit_y_uvec);

  if((!(index_NA.n_elem==0)) && ((index_full.n_elem>15 && s=="lda")|| (index_full.n_elem>posit_x.n_elem && s!="lda"))){

  //dividing data to NA and full

  arma::mat X_full = x(index_full,posit_x);
  arma::mat X_NA = x(index_NA,posit_x);
  arma::colvec Y_full = x(index_full,posit_y_uvec);
  arma::colvec w_full = w.elem(index_full);

  pfuncw f = funMapw[s];
  pred  = (*f)(Y_full,X_full,w_full,X_NA,times);


  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;

  return Y;
}

//eval vifs


arma::mat cov2cor(arma::mat & X){
  arma::mat corr = arma::diagmat(1/arma::sqrt(X.diag())) * X * arma::diagmat(1/arma::sqrt(X.diag()));
  return corr;
}


// [[Rcpp::export]]
arma::vec VIF_(arma::mat &x,int posit_y,arma::uvec posit_x,arma::uvec posit_x_var,bool correct){

  arma::mat x_cols = x.cols(posit_x - 1);

  arma::colvec y_cols = x.col(posit_y - 1);

  arma::uvec Nvar_x_o = arma::unique(posit_x_var);
  int Nvar_o = Nvar_x_o.n_elem;

  arma::vec vifs(Nvar_o,arma::fill::none);


  arma::uvec full_rows = get_index_full_R(x,posit_y - 1,posit_x - 1);

  arma::mat full_x = x_cols.rows(full_rows);

  arma::mat full_y = y_cols.rows(full_rows);


  arma::uvec vols = arma::find(arma::var(full_x)>0);

  int vols_n = vols.n_elem;


  arma::mat x_vols = full_x.cols(vols);


  if(vols_n<Nvar_o){Rcpp::warning("There is at least a one zero variance variable");return vifs;}


  arma::uvec posit_x_v = posit_x.elem(vols);

  arma::uvec posit_x_var_v = posit_x_var.elem(vols);

  arma::mat col_means = arma::mean(x_vols,0);

  x_vols.each_row() -= col_means;


  int Ncol = posit_x_v.n_elem;

  arma::uvec Nvar_x = arma::unique(posit_x_var_v);
  int Nvar = Nvar_x.n_elem;


  arma::mat xtx = x_vols.t()*x_vols;

  arma::mat XXinv = arma::inv_sympd(xtx);


  if(Nvar<Ncol){ XXinv = cov2cor(XXinv) ;}


  double det_XXinv = arma::det(XXinv);


  arma::vec dfs(Nvar);

  arma::uvec pp = arma::linspace<arma::uvec>(0, Ncol-1, Ncol);

  for(int i=0;i<(Nvar);i++){

    int a = Nvar_x(i);

    arma::mat XXinv_small = XXinv;

    arma::uvec ii = pp.elem(arma::find(posit_x_var_v==a));

    int ii_len = ii.n_elem;

    int ii_first = ii(0);

    int ii_last = ii(ii_len-1);

    dfs(i) = ii_len;

    XXinv_small.shed_rows(ii_first,ii_last);

    XXinv_small.shed_cols(ii_first,ii_last);

    double vif = arma::as_scalar(arma::det(XXinv(ii,ii)) * arma::det(XXinv_small)/det_XXinv);

    vifs(i) = vif;

  }

  if(correct){
    for(unsigned int i=0;i<(vifs.n_elem);i++){
      vifs(i) = pow(vifs(i),(1/(2*dfs(i))));
    }
  }

  return vifs;

};



// [[Rcpp::export]]
arma::colvec fill_NA_N_(arma::mat &x, std::string model, int posit_y,arma::uvec posit_x,arma::colvec w,int times=10){


  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w.is_empty() || (model.compare("lda") == 0)){
    pred_avg = impute_raw_R(x,model,posit_y,posit_x,times);
  } else{
    pred_avg = imputeW_R(x,model,posit_y,posit_x,w,times);
  }

  //index
  return pred_avg;
}


// [[Rcpp::export]]
arma::colvec fill_NA_(arma::mat &x,std::string model, int posit_y,arma::uvec posit_x,arma::colvec w){

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w.is_empty() || (model.compare("lda") == 0)){
    pred_avg = impute_raw_R(x,model,posit_y,posit_x,1);
  } else {
    pred_avg = imputeW_R(x,model,posit_y,posit_x,w,1);
  }

  //index
  return pred_avg;
}

// // [[Rcpp::export]]
// arma::colvec PMM(arma::mat &x, int posit_y,arma::uvec posit_x,arma::colvec w, int k){
//
//   posit_x =  posit_x - 1;
//   posit_y = posit_y - 1;
//
//   arma::uvec posit_y_uvec(1);
//   posit_y_uvec(0) = posit_y;
//
//   arma::uvec index_full = get_index_full_R(x,posit_y, posit_x);
//   arma::uvec index_NA = get_index_NA_R(x,posit_y, posit_x);
//
//   arma::colvec pred = x(index_NA,posit_y_uvec);
//
//   if((!(index_NA.n_elem==0)) && ((index_full.n_elem>15 && s=="lda")|| (index_full.n_elem>posit_x.n_elem && s!="lda"))){
//
//     arma::mat X_full = x(index_full,posit_x);
//     arma::mat X_NA = x(index_NA,posit_x);
//     arma::colvec Y_full = x(index_full,posit_y_uvec);
//
//   int N = X.n_rows; int C = X.n_cols; int N_NA = X1.n_rows;
//   double df = N-C ;
//
//   arma::colvec coef(C);
//
//   if(w.is_empty()){
//
//     coef = arma::solve(X.t()*X, X.t()*y);
//
//   } else {
//
//     arma::colvec wq=sqrt(w);
//     arma::colvec y2 = wq%y;
//     arma::mat X2(N,C);
//
//     for(int h=0;h<C;h++){
//       X2.col(h)=wq%X.col(h);
//     }
//
//     coef = arma::solve(X2.t()*X2, X2.t()*y2);
//
//
//   }
//
//   arma::colvec res = y - X*coef;
//
//   //arma::mat XX_inv = arma::inv(arma::trans(X)*X) ;
//
//   double res2 = arma::as_scalar(arma::trans(res)*res);
//
//   double chi2 = Rcpp::as<double>(Rcpp::rchisq(1, df));
//
//   double sigma_b = sqrt(res2/chi2);
//
//   arma::vec noise2(N_NA);
//   noise2.randn();
//
//   arma::colvec miss(N_NA,arma::fill::zeros);
//
//   miss = (X1 * coef + noise2 * sigma_b);
//
//
//
//
//   //index
//   return pred_avg;
// }
