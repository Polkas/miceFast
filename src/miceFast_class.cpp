//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>
#include <map>
#include "miceFast.h"

//miceFast Constructors


miceFast::miceFast(arma::mat _x):
  x(_x){};

miceFast::miceFast(arma::mat _x,arma::uvec _g,bool _s):
  x(_x),
  g(_g),
  sorted(_s){};

miceFast::miceFast(arma::mat _x,arma::colvec _w):
  x(_x),
  w(_w){};

miceFast::miceFast(arma::mat _x,arma::uvec _g,bool _s,arma::colvec _w):
  x(_x),
  g(_g),
  sorted(_s),
  w(_w){};


//functions for getting indexes

arma::uvec miceFast::get_index_full(int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  arma::uvec index_Full = arma::find((full_y + full_x) == 2);

  return index_Full;
}

arma::uvec miceFast::get_index_NA(int posit_y, arma::uvec posit_x){

  arma::colvec y_col = x.col(posit_y);

  arma::mat x_cols = x.cols(posit_x);

  arma::uvec full_y = complete_cases_vec(y_col);

  arma::uvec full_x = complete_cases_mat(x_cols);

  int N_rows = x.n_rows;

  arma::uvec index(N_rows,arma::fill::zeros);

  for(int i=0;i<N_rows;i++){

    if( (full_y(i)==0) && (full_x(i)==1) ){index(i) = 1;}

  }

  arma::uvec index_NA = arma::find(index==1);


  return index_NA;
}


arma::uvec miceFast::get_index_full_R(int posit_y, arma::uvec posit_x){

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::uvec index_f = miceFast::get_index_full(posit_y,posit_x) + 1L;

  arma::uvec index_full = index_f;

  return index_full;

}

arma::uvec miceFast::get_index_NA_R(int posit_y, arma::uvec posit_x){

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::uvec index_n = miceFast::get_index_NA(posit_y,posit_x) + 1L;

  arma::uvec index_NA = index_n;

  return index_NA;
}

void miceFast::sortData_byg(){

  arma::uvec order = arma::stable_sort_index(g);
  x = x.rows(order);
  g = g.elem(order);
  if(!w.is_empty()){w = w.elem(order);}

}


//function for sorting data by a grouping variable

bool miceFast::is_vars_updated(){

  return updated;

}

//function for printing recommended prediction models

std::string miceFast::get_models(int posit_y){

  arma::colvec Y_raw = x.col(posit_y - 1);

  arma::uvec index_full =  arma::find_finite(Y_raw);

  if(!Y_raw.has_nan()){return "no NA values for the dependent variable";}

  arma::colvec  Y = Y_raw.rows(index_full);

  arma::vec un = arma::unique(Y);


  int un_n =  un.n_elem;
  std::string type;

  if(un_n == 2){
    type = "recommended lda or (lm_pred,lm_bayes,lm_noise - round results)";
  }
  else if(un_n > 2 && un_n <= 15){
    type = "recommended lda or (lm_pred,lm_bayes,lm_noise - remember to round results if needed)";
  }
  else if(un_n > 15){
    type = "lm_pred or lm_bayes or lm_noise";
  }
  else {type = "one unique value";}

  return type;
}

std::string miceFast::get_model(int posit_y){

  arma::colvec Y_raw = x.col(posit_y - 1);

  arma::uvec index_full =  arma::find_finite(Y_raw);

  if(!Y_raw.has_nan()){return "no NA values for the dependent variable";}

  arma::colvec  Y = Y_raw.rows(index_full);

  arma::vec un = arma::unique(Y);

  int un_n =  un.n_elem;
  std::string type;

  if(un_n == 2){
    type = "lda";
  }
  else if(un_n > 2 && un_n <= 15){
    type = "lda";
  }
  else if(un_n > 15){
    type = "lm_pred";
  }
  else {type = "one unique value";}

  return type;
}

// map - implementing functions

typedef arma::colvec (*pfunc)(arma::colvec&,arma::mat&,arma::mat&);
std::map<std::string, pfunc> funMap = {{"lda",fastLda},
{"lm_pred",fastLm_pred},
{"lm_noise",fastLm_noise},
{"lm_bayes",fastLm_bayes}};

//Impute


// arma::mat miceFast::impute_auto(){
//
//   unsigned int nrows = x.n_rows;
//   unsigned int ncols = x.n_cols;
//
//   arma::mat new_x = arma::mat(nrows,ncols,arma::fill::none);
//
//   for(unsigned int i=0; i < (ncols - 1)  ;i++){
//
//     std::string model_type = miceFast::get_model(i + 1);
//
//     arma::colvec pred(nrows);
//
//     if(model_type == "one unique value" || model_type == "no NA values for the dependent variable"){
//
//       pred = x.col(i);
//
//     } else {
//       arma::uvec cols = arma::regspace<arma::uvec>(0,(ncols - 1));
//
//       cols.shed_row(i);
//
//       Rcpp::List preds =  miceFast::impute(model_type ,i + 1, cols + 1 , true);
//
//       Rcpp::NumericVector p = preds["imputations"];
//
//       arma::colvec pred(p.begin(),p.size(), false);
//
//       }
//
//     new_x.col(i) = pred;
//   }
//
//   return new_x;
// }


Rcpp::List miceFast::impute_force(std::string s, int posit_y,arma::uvec posit_x){

  updated = true;

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred;

  if(w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputeby(s,posit_y,posit_x);
  }
  else if(w.is_empty() && g.is_empty())
  {
    pred = miceFast::impute_raw(s,posit_y,posit_x);
  }
  else if(!w.is_empty() && g.is_empty())
  {
    pred = miceFast::imputeW(s,posit_y,posit_x);
  }
  else if(!w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputebyW(s,posit_y,posit_x);
  }

  x.col(posit_y) = pred;

  arma::uvec index_NA = miceFast::get_index_NA_R(posit_y+1,posit_x+1);
  arma::uvec index_full = miceFast::get_index_NA_R(posit_y+1,posit_x+1);

  return Rcpp::List::create(Rcpp::Named("imputations") = pred,
                            Rcpp::Named("index_NA") = index_NA,
                            Rcpp::Named("index_full") = index_full);

}


Rcpp::List miceFast::impute(std::string s, int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred;

  if(w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputeby(s,posit_y,posit_x);
  }
  else if(w.is_empty() && g.is_empty())
  {
    pred = miceFast::impute_raw(s,posit_y,posit_x);
  }
  else if(!w.is_empty() && g.is_empty())
  {
    pred = miceFast::imputeW(s,posit_y,posit_x);
  }
  else if(!w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputebyW(s,posit_y,posit_x);
  }

  arma::uvec index_NA = miceFast::get_index_NA_R(posit_y+1,posit_x+1);
  arma::uvec index_full = miceFast::get_index_NA_R(posit_y+1,posit_x+1);

  return Rcpp::List::create(Rcpp::Named("imputations") = pred,
                            Rcpp::Named("index_NA") = index_NA,
                            Rcpp::Named("index_full") = index_full);

}

arma::colvec miceFast::impute_raw(std::string s, int posit_y,arma::uvec posit_x){

  arma::uvec index_full = miceFast::get_index_full(posit_y, posit_x);


  arma::uvec index_NA = miceFast::get_index_NA(posit_y, posit_x);

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.rows(index_full);

  pfunc f = funMap[s];
  arma::colvec pred = (*f)(Y_full,X_full,X_NA);

  Y.rows(index_NA) = pred;

  return Y;
}


//Impute with grouping

arma::colvec miceFast::imputeby(std::string s, int posit_y,arma::uvec posit_x){

  if(sorted == false){
    sortData_byg();
  }

  //index

  arma::uvec index_full = get_index_full(posit_y,posit_x);
  arma::uvec index_NA = get_index_NA(posit_y,posit_x);

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  //if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(g.has_nan()){Rcpp::stop("There is NA values for the grouping variable");}

  //grouping variable

  arma::uvec un = arma::unique(g);

  unsigned int group = un.n_elem;

  //quantitative model

  pfunc fun = funMap[s];

  //dividing data to NA and full

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.rows(index_full);

  arma::uvec g_full = g.elem(index_full);
  arma::uvec g_NA = g.elem(index_NA);

  //predictions container

  arma::colvec pred_all(index_NA.n_elem);

  //iter

  arma::uvec his_full = histFast(g_full);
  arma::uvec ends_full = arma::cumsum(his_full)-1;
  arma::uvec starts_full = arma::shift(ends_full,1)+1;
  starts_full(0) = 0;

  arma::uvec his_NA = histFast(g_NA);
  arma::uvec ends_NA = arma::cumsum(his_NA)-1;
  arma::uvec starts_NA = arma::shift(ends_NA,1)+1;
  starts_NA(0) = 0;

  unsigned int start = 0;
  unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    arma::mat X_full_0 = X_full.rows(starts_full(a),ends_full(a)) ;  // copying. It would be better to improve it to use reference to certain rows
    arma::mat X_NA_0 = X_NA.rows(starts_NA(a),ends_NA(a));  // copying. It would be better to improve it to use reference to certain rows
    arma::colvec Y_full_0 =  Y_full.subvec(starts_full(a),ends_full(a));  // copying. It would be better to improve it to use reference to certain rows

    arma::colvec pred =  (*fun)(Y_full_0,X_full_0,X_NA_0);

    end = start + pred.n_elem-1;

    pred_all(arma::span(start,end)) = pred;

    start = end+1;

  }

  Y.rows(index_NA) = pred_all;

  return Y;
}

//WEIGHTED

typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&);
std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
{"lm_noise",fastLm_weighted_noise},
{"lm_bayes",fastLm_weighted_bayes}};

arma::colvec miceFast::imputeW(std::string s,int posit_y,arma::uvec posit_x){

  //index

  arma::uvec index_full = get_index_full(posit_y,posit_x);
  arma::uvec index_NA = get_index_NA(posit_y,posit_x);

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights");}
  if(arma::any(w<0)){Rcpp::stop("There are ngative values for the weights variable");}

  //dividing data to NA and full

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.elem(index_full);
  arma::colvec w_full = w.elem(index_full);

  pfuncw f = funMapw[s];
  arma::colvec pred = (*f)(Y_full,X_full,w_full,X_NA);

  Y.rows(index_NA) = pred;

  return Y;
}

//Impute with grouping

arma::colvec miceFast::imputebyW(std::string s,int posit_y,arma::uvec posit_x){

  if(sorted == false){
    sortData_byg();
  }

  //index

  arma::uvec index_full = get_index_full(posit_y,posit_x);
  arma::uvec index_NA = get_index_NA(posit_y,posit_x);

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights variable");}
  if(arma::any(w<0)){Rcpp::stop("There are ngative values for the weights variable");}
  if(g.has_nan()){Rcpp::stop("There is NA values for the grouping variable");}

  //grouping variable

  arma::uvec un = arma::unique(g);

  unsigned int group = un.n_elem;

  //quantitative model

  pfuncw fun = funMapw[s];

  //dividing data to NA and full

  arma::mat X_full = X.rows(index_full);
  arma::mat X_NA = X.rows(index_NA);
  arma::colvec Y_full = Y.rows(index_full);
  arma::colvec w_full = w.rows(index_full);


  arma::uvec g_full = g.elem(index_full);
  arma::uvec g_NA = g.elem(index_NA);

  //predictions container

  arma::colvec pred_all(index_NA.n_elem);

  // start end

  //iter

  arma::uvec his_full = histFast(g_full);
  arma::uvec ends_full = arma::cumsum(his_full)-1;
  arma::uvec starts_full = arma::shift(ends_full,1)+1;
  starts_full(0) = 0;

  arma::uvec his_NA = histFast(g_NA);
  arma::uvec ends_NA = arma::cumsum(his_NA)-1;
  arma::uvec starts_NA = arma::shift(ends_NA,1)+1;
  starts_NA(0) = 0;


  unsigned int start = 0;
  unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    arma::mat X_full_0 = X_full.rows(starts_full(a),ends_full(a)) ;  // copying. It would be better to improve it to use reference to certain rows
    arma::mat X_NA_0 = X_NA.rows(starts_NA(a),ends_NA(a));  // copying. It would be better to improve it to use reference to certain rows
    arma::colvec Y_full_0 =  Y_full.subvec(starts_full(a),ends_full(a));  // copying. It would be better to improve it to use reference to certain rows
    arma::colvec w_full_0 =  w_full.subvec(starts_full(a),ends_full(a)); // copying. It would be better to improve it to use reference to certain rows

    arma::colvec pred =  (*fun)(Y_full_0,X_full_0,w_full_0,X_NA_0);

    end = start + pred.n_elem-1;

    pred_all(arma::span(start,end)) = pred;

    start = end+1;

  }

  Y.rows(index_NA) = pred_all;

  return Y;
}

RCPP_MODULE(miceFast){
  using namespace Rcpp ;

  class_<miceFast>("miceFast")
    .constructor<arma::mat>()
    .constructor<arma::mat,arma::uvec,bool>()
    .constructor<arma::mat,arma::colvec>()
    .constructor<arma::mat,arma::uvec,bool,arma::colvec>()
    .method("impute", &miceFast::impute)
    .method("impute_force", &miceFast::impute_force)
    .method("get_models", &miceFast::get_models)
    .method("get_model", &miceFast::get_model)
    .method("is_vars_updated", &miceFast::is_vars_updated)
    .method("get_index_NA_R", &miceFast::get_index_NA_R)
    .method("get_index_full_R", &miceFast::get_index_full_R)
  ;}
