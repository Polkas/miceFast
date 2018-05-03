#include <RcppArmadillo.h>
#include "miceFast.h"

//functions for getting indexes

arma::uvec get_index_full(arma::mat &x,int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index_Full = arma::find((full_y + full_x) == 2);

  return index_Full;
}

arma::uvec get_index_NA(arma::mat &x,int posit_y, arma::uvec posit_x){

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

arma::colvec impute_raw(arma::mat &x,std::string s, int posit_y,arma::uvec posit_x,int times){

  typedef arma::colvec (*pfunc)(arma::colvec&,arma::mat&,arma::mat&,int);
  std::map<std::string, pfunc> funMap = {
    {"lda",fastLda},
    {"lm_pred",fastLm_pred},
    {"lm_noise",fastLm_noise},
    {"lm_bayes",fastLm_bayes}};

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  arma::uvec index_full = get_index_full(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA(x,posit_y, posit_x);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.rows(index_full);

  pfunc f = funMap[s];
  arma::colvec pred = (*f)(Y_full,X_full,X_NA,times);

  Y.rows(index_NA) = pred;

  return Y;
}


arma::colvec imputeW(arma::mat &x,std::string s,int posit_y,arma::uvec posit_x,arma::colvec w,int times){

  typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&,int);
  std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
  {"lm_noise",fastLm_weighted_noise},
  {"lm_bayes",fastLm_weighted_bayes}};

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  arma::uvec index_full = get_index_full(x,posit_y, posit_x);
  arma::uvec index_NA = get_index_NA(x,posit_y, posit_x);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights");}
  if(arma::any(w<0)){Rcpp::stop("There are negative values for the weights variable");}

  //dividing data to NA and full

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.elem(index_full);
  arma::colvec w_full = w.elem(index_full);

  pfuncw f = funMapw[s];
  arma::colvec pred = (*f)(Y_full,X_full,w_full,X_NA,times);

  Y.rows(index_NA) = pred;

  return Y;
}

//eval vif

//' \code{vif}
//'
//' @param x a numeric matrix - a numeric matrix with variables
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//'
//' @return load a numeric vector with VIF for all variables matching with a provided posit_x
//'
//' @seealso \code{\link{impute_N}}
//'
//' @examples
//' \dontrun{
//'    impute(filename)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec vifs(arma::mat &x,int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  arma::mat x_cols = x.cols(posit_x - 1);

  arma::uvec full_rows = get_index_full(x,posit_y - 1,posit_x - 1);

  arma::mat full_x = x_cols.rows(full_rows);

  arma::mat col_means = arma::mean(full_x,0);

  full_x.each_row() -= col_means;

  arma::mat XXinv = arma::inv(full_x.t()*full_x);

  int Ncol = full_x.n_cols;

  arma::vec vifs(Ncol);

  double det_XXinv = arma::det(XXinv);

  for(int i=0;i<Ncol;i++){

    arma::mat XXinv_small = XXinv;

    XXinv_small.shed_row(i);

    XXinv_small.shed_col(i);

    vifs(i)= XXinv(i,i) * arma::det(XXinv_small)/det_XXinv;

  }

  return vifs;
};


//' \code{impute_N}
//'
//' @param dt data.table - a data.table with variables
//' @param model a character - a posibble options ("lm_bayes","lm_noise")
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//' @param w  a numeric vector - a weighting variable - only positive values
//' @param times an integer - a number of multiple imputations
//'
//' @return load data in a data.table format
//'
//' @seealso \code{\link{impute}}
//'
//' @examples
//' \dontrun{
//'    fars_read(filename)
//' }
//'
//' @export
//'
// [[Rcpp::export]]
arma::colvec impute_N(arma::mat &x, std::string s, int posit_y,arma::uvec posit_x,arma::colvec w = 0,int times){

  if( !(s.compare("lm_bayes") == 0) && !(s.compare("lm_noise") == 0)){Rcpp::stop("Works only for `lm_bayes` and `lm_noise` models");}
  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w.n_elem<2){
   pred_avg = impute_raw(x,s,posit_y,posit_x,times);
  } else {
    pred_avg = imputeW(x,s,posit_y,posit_x,w,times);
  }

  //index

  return pred_avg;

}

//' \code{impute}
//'
//' @param dt data.table - a data.table with variables
//' @param model a character - a posibble options ("lda","lm_pred","lm_bayes","lm_noise")
//' @param posit_y an integer - a position of dependent variable
//' @param posit_x an integer vector - positions of independent variables
//' @param w  a numeric vector - a weighting variable - only positive values
//'
//' @return load data in a data.table format
//'
//' @seealso \code{\link{impute}}
//'
//' @examples
//' \dontrun{
//'    fars_read(filename)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::colvec impute(arma::mat &x,std::string s, int posit_y,arma::uvec posit_x,arma::colvec w = 0){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg;

  if(w.n_elem<2){
    pred_avg = impute_raw(x,s,posit_y,posit_x,1);
  } else {
    pred_avg = imputeW(x,s,posit_y,posit_x,w,1);
  }

  //index

  return pred_avg;


}

