//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>
#include <map>
#include "miceFast.h"

//miceFast Constructors

miceFast::miceFast(){};

miceFast::~miceFast(){

  x.clear();
  g.clear();
  w.clear();
  index.clear();
  sorted=false;
  updated.clear();

};

void miceFast::set_data(arma::mat & _x){
  x = arma::mat(_x.begin(),_x.n_rows,_x.n_cols,false);
  index = arma::regspace<arma::uvec>(0,x.n_rows - 1) + 1;
}

void miceFast::set_g(arma::uvec & _g){
  g =_g;
  sorted = g.is_sorted();

}

void miceFast::set_w(arma::colvec & _w){
  w = arma::colvec(_w.begin(),_w.n_rows,false);;
}

bool miceFast::is_sorted_byg(){
  if(g.is_empty()){Rcpp::stop("There is no grouping variable provided");}
  return sorted;
};

arma::mat miceFast::get_data(){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  return x;
};

arma::colvec miceFast::get_w(){
  if(w.is_empty()){Rcpp::stop("There is no weighting variable provided");}
  return w;
};
arma::uvec miceFast::get_g(){
  if(g.is_empty()){Rcpp::stop("There is no grouping variable provided");}
  return g;
};

arma::uvec miceFast::get_index(){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  return index;
};

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


void miceFast::sortData_byg(){

  if(g.is_empty()){Rcpp::stop("There is no grouping variable provided");}
  Rcpp::warning("\n Data was sorted by the grouping variable - use `get_index()` to retrieve an order");
  arma::uvec order = arma::stable_sort_index(g);
  x = x.rows(order);
  g = g.elem(order);
  index = index.elem(order);
  if(!w.is_empty()){w = w.elem(order);}
  sorted = true;

}

//function for sorting data by a grouping variable

std::vector<int> miceFast::which_updated(){

  return updated;

}

//function for printing recommended prediction models

std::string miceFast::get_models(int posit_y){

  if(x.is_empty()){Rcpp::stop("at least set the data");}

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

  if(x.is_empty()){Rcpp::stop("at least set the data");}

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


void miceFast::update_var(int posit_y,arma::vec impute){

  if(x.is_empty()){Rcpp::stop("at least set the data");}
  if( x.n_rows != impute.n_elem ){Rcpp::stop("wrong number of observations");}

  posit_y = posit_y - 1;

  x.col(posit_y) = impute;

  updated.push_back(posit_y + 1);

}


Rcpp::List miceFast::impute(std::string s, int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  Rcpp::List pred =  option_impute(s,posit_y,posit_x);

  return Rcpp::List::create(Rcpp::Named("imputations") = pred["imputations"],
                            Rcpp::Named("index_imputed") = pred["index_NA"],
                            Rcpp::Named("index_full") = pred["index_full"]);

}

Rcpp::List miceFast::option_impute(std::string s,int posit_y,arma::uvec posit_x){

  Rcpp::List pred;

  if(w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputeby(s,posit_y,posit_x);
  }
  else if(w.is_empty() && g.is_empty())
  {
    pred = miceFast::impute_raw(s,posit_y,posit_x);
  }
  else if( !w.is_empty() && g.is_empty())
  {
    if(s=="lda"){pred = miceFast::impute_raw(s,posit_y,posit_x);} else {pred = miceFast::imputeW(s,posit_y,posit_x);}
  }
  else if(!w.is_empty() && !g.is_empty())
  {
    if(s=="lda"){pred = miceFast::imputeby(s,posit_y,posit_x);} else {pred = miceFast::imputebyW(s,posit_y,posit_x);}
  }
  return pred;

}

Rcpp::List miceFast::impute_raw(std::string s, int posit_y,arma::uvec posit_x){

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

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);
  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);
  index_full_return.elem(index_full).fill(1);

  return Rcpp::List::create(Rcpp::Named("imputations") = Y,
                            Rcpp::Named("index_NA") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);
}


//Impute with grouping

Rcpp::List miceFast::imputeby(std::string s, int posit_y,arma::uvec posit_x){

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

  arma::uvec his_full = arma::hist(g_full,un);
  arma::uvec ends_full = arma::cumsum(his_full)-1;
  arma::uvec starts_full = arma::shift(ends_full,1) + 1;
  starts_full(0) = 0;

  arma::uvec his_NA = arma::hist(g_NA,un);
  arma::uvec ends_NA = arma::cumsum(his_NA)-1;
  arma::uvec starts_NA = arma::shift(ends_NA,1) + 1;
  starts_NA(0) = 0;

  unsigned int start = 0;
  unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    int ss_NA = starts_NA(a);
    int ee_NA = ends_NA(a);
    int ss_full = starts_full(a);
    int ee_full = ends_full(a);

    if((ss_NA < ee_NA) && (ss_full <= ee_full)){

      arma::mat X_full_0 = X_full.rows(ss_full,ee_full) ;  // copying. It would be better to improve it to use reference
      arma::mat X_NA_0 = X_NA.rows(ss_NA,ee_NA);  // copying. It would be better to improve it to use reference
      arma::colvec Y_full_0 =  Y_full.rows(ss_full,ee_full);  // copying. It would be better to improve it to use reference

      arma::colvec pred =  (*fun)(Y_full_0,X_full_0,X_NA_0);

      end = start + pred.n_elem - 1;

      pred_all.rows(start,end) = pred;

      start = end + 1;

    }
  }

  Y.rows(index_NA) = pred_all;

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);
  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);
  index_full_return.elem(index_full).fill(1);

  return Rcpp::List::create(Rcpp::Named("imputations") = Y,
                            Rcpp::Named("index_NA") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);
}

//WEIGHTED

typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&);
std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
{"lm_noise",fastLm_weighted_noise},
{"lm_bayes",fastLm_weighted_bayes}};

Rcpp::List miceFast::imputeW(std::string s,int posit_y,arma::uvec posit_x){

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

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);
  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);
  index_full_return.elem(index_full).fill(1);

  return Rcpp::List::create(Rcpp::Named("imputations") = Y,
                            Rcpp::Named("index_NA") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);
}

//Impute with grouping

Rcpp::List miceFast::imputebyW(std::string s,int posit_y,arma::uvec posit_x){

  if(!sorted){
    sortData_byg();
  }

  //index
  arma::uvec index_full = get_index_full(posit_y,posit_x);
  arma::uvec index_NA = get_index_NA(posit_y,posit_x);

  arma::mat X = x.cols(posit_x);
  arma::colvec Y = x.col(posit_y);

  if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights variable");}
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

  arma::uvec his_full = arma::hist(g_full,un);
  arma::uvec ends_full = arma::cumsum(his_full)-1;
  arma::uvec starts_full = arma::shift(ends_full,1) + 1;
  starts_full(0) = 0;

  arma::uvec his_NA = arma::hist(g_NA,un);
  arma::uvec ends_NA = arma::cumsum(his_NA) - 1;
  arma::uvec starts_NA = arma::shift(ends_NA,1) + 1;
  starts_NA(0) = 0;


  unsigned int start = 0;
  unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    int ss_NA = starts_NA(a);
    int ee_NA = ends_NA(a);
    int ss_full = starts_full(a);
    int ee_full = ends_full(a);

    if((ss_NA < ee_NA) && (ss_full <= ee_full)){

    arma::mat X_full_0 = X_full.rows(ss_full,ee_full) ;  // copying. It would be better to improve it to use reference
    arma::mat X_NA_0 = X_NA.rows(ss_NA,ee_NA);  // copying. It would be better to improve it to use reference
    arma::colvec Y_full_0 =  Y_full.rows(ss_full,ee_full);  // copying. It would be better to improve it to use reference
    arma::colvec w_full_0 =  w_full.rows(ss_full,ee_full); // copying. It would be better to improve it to use reference

    arma::colvec pred =  (*fun)(Y_full_0,X_full_0,w_full_0,X_NA_0);

    end = start + pred.n_elem - 1;

    pred_all.rows(start,end) = pred;

    start = end + 1;

    }

  }

  Y.rows(index_NA) = pred_all;

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);
  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);
  index_full_return.elem(index_full).fill(1);

  return Rcpp::List::create(Rcpp::Named("imputations") = Y,
                            Rcpp::Named("index_NA") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);
}

RCPP_MODULE(miceFast){
  using namespace Rcpp ;

  class_<miceFast>("miceFast")
    .default_constructor()
    //.field("data",&miceFast::x)
    //.field("g",&miceFast::g)
    //.field("w",&miceFast::w)
    //.field("updated",&miceFast::updated)
    //.field("sorted",&miceFast::sorted)
    .method("get_data", &miceFast::get_data)
    .method("get_w", &miceFast::get_w)
    .method("get_g", &miceFast::get_g)
    .method("get_index", &miceFast::get_index)
    .method("is_sorted_byg", &miceFast::is_sorted_byg)
    .method("which_updated", &miceFast::which_updated)

    .method("impute", &miceFast::impute)
    .method("update_var", &miceFast::update_var)
    .method("get_models", &miceFast::get_models)
    .method("get_model", &miceFast::get_model)
    .method("set_data", &miceFast::set_data)
    .method("set_g", &miceFast::set_g)
    .method("set_w", &miceFast::set_w)
    .method("sort_byg",&miceFast::sortData_byg)

  ;}
