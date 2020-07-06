//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>
#include <map>
#include "miceFast.h"

//miceFast Constructors

miceFast::miceFast(){}

miceFast::~miceFast(){

//x.clear();
//g.clear();
//w.clear();
//index.clear();
//sorted=false;
//updated.clear();

}

void miceFast::set_data(arma::mat & _x){
  x = arma::mat(_x.begin(),_x.n_rows,_x.n_cols,false,false);
  N_rows = x.n_rows;
  N_cols = x.n_cols;
  index = arma::regspace<arma::uvec>(0,N_rows - 1) + 1;
}

void miceFast::set_data_sparse(arma::sp_mat & _x){
  x_sp = _x;
  N_rows = x_sp.n_rows;
  N_cols = x_sp.n_cols;
  index = arma::regspace<arma::uvec>(0,N_rows - 1) + 1;
}

void miceFast::set_g(arma::colvec & _g){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  unsigned int n_elems = _g.n_rows;
  if(N_rows!=n_elems){Rcpp::stop("Wrong number of elements");}
  g = arma::colvec(_g.begin(),n_elems,false,false);  //g = arma::uvec(_g.begin(),_g.n_elem,false,false); - not work for uvec - unit not equal to int
  sorted = g.is_sorted();
}

void miceFast::set_w(arma::colvec & _w){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  unsigned int n_elems = _w.n_rows;
  if(N_rows!=n_elems){Rcpp::stop("Wrong number of elements");}
  w = arma::colvec(_w.begin(),n_elems,false,false);
}

void miceFast::set_ridge(double _ridge){
  ridge = _ridge;
}
bool miceFast::is_sorted_byg(){
  if(g.is_empty()){Rcpp::stop("There is no grouping variable provided");}
  return sorted;
}

arma::mat miceFast::get_data(){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  return this -> x;
}

arma::colvec miceFast::get_w(){
  if(w.is_empty()){Rcpp::stop("There is no weighting variable provided");}
  return w;
}
arma::colvec miceFast::get_g(){
  if(g.is_empty()){Rcpp::stop("There is no grouping variable provided");}
  return g;
}

double miceFast::get_ridge(){
  return ridge;
}

arma::uvec miceFast::get_index(){
  if(x.is_empty()){Rcpp::stop("There is no data provided");}
  return index;
}

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

  arma::uvec index(N_rows,arma::fill::zeros);

  for(unsigned int i=0;i<N_rows;i++){

    if( (full_y(i)==0) && (full_x(i)==1) ){index(i) = 1;}

  }

  arma::uvec index_NA = arma::find(index==1);

  return index_NA;
}


void miceFast::sortData_byg(){

  if(g.is_empty()){Rcpp::stop("There is no a grouping variable provided");}
  if(!sorted){

  Rcpp::warning("\n Data was sorted by the grouping variable - use `get_index()` to retrieve an order");

  arma::uvec order = arma::stable_sort_index(g);

  arma::mat x_temp = x.rows(order);
  x = x_temp;
  x_temp.clear();

  g.col(0) = g.rows(order);

  index = index.elem(order);

  if(!w.is_empty()){w.col(0) =  w.rows(order);}

  sorted = true;

  }

}

//function for sorting data by a grouping variable

std::vector<int> miceFast::which_updated(){

  return updated;

}

//eval vif
arma::vec miceFast::vifs(int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  arma::mat x_cols = x.cols(posit_x - 1);

  arma::uvec full_rows = get_index_full(posit_y - 1,posit_x - 1);

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


void miceFast::update_var(int posit_y,arma::vec impute){

  if(x.is_empty()){Rcpp::stop("at least set the data");}
  if( N_rows != impute.n_elem ){Rcpp::stop("wrong number of observations");}

  posit_y = posit_y - 1;

  x.col(posit_y) = impute;

  updated.push_back(posit_y + 1);

}

Rcpp::List miceFast::impute_N(std::string s, int posit_y,arma::uvec posit_x,int k){

  if( !(s.compare("lm_bayes") == 0) && !(s.compare("lm_noise") == 0) && !(s.compare("pmm") == 0)){Rcpp::stop("Works only for `lm_bayes`, `lm_noise` and `pmm` models");}
  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  arma::colvec pred_avg =  option_impute_multiple(s,posit_y,posit_x,k);

  //index

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);

  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);

  index_full_return.elem(index_full).fill(1);

  return Rcpp::List::create(Rcpp::Named("imputations") = pred_avg,
                            Rcpp::Named("index_imputed") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);

}

Rcpp::List miceFast::impute(std::string s, int posit_y,arma::uvec posit_x){

  if(!different_y_and_x(posit_y,posit_x)){Rcpp::stop("the same variable is dependent and indepentent");}
  if(!different_x(posit_x)){Rcpp::stop("the same variables repeated few times as independent");}
  if(x.is_empty()){Rcpp::stop("at least set the data");}

  posit_x =  posit_x - 1;
  posit_y = posit_y - 1;

  const int k = 1;

  arma::colvec pred =  option_impute_multiple(s,posit_y,posit_x,k);

  //index

  arma::uvec index_NA_return(x.n_rows,arma::fill::zeros);

  index_NA_return.elem(index_NA).fill(1);

  arma::uvec index_full_return(x.n_rows,arma::fill::zeros);

  index_full_return.elem(index_full).fill(1);


  return Rcpp::List::create(Rcpp::Named("imputations") = pred,
                            Rcpp::Named("index_imputed") = index_NA_return,
                            Rcpp::Named("index_full") = index_full_return);

}

// map - implementing functions

  typedef arma::colvec (*pfunc)(arma::colvec&,arma::mat&,arma::mat&,int, double);
    std::map<std::string, pfunc> funMap = {
    {"lda",fastLda},
    {"lm_pred",fastLm_pred},
    {"lm_noise",fastLm_noise},
    {"lm_bayes",fastLm_bayes},
    {"pmm",pmm_neibo}};



arma::colvec miceFast::option_impute_multiple(std::string s,int posit_y,arma::uvec posit_x,int k){

  arma::colvec pred;

  if(w.is_empty() && !g.is_empty())
  {
    pred = miceFast::imputeby(s,posit_y,posit_x,k);
  }
  else if(w.is_empty() && g.is_empty())
  {
    pred = miceFast::impute_raw(s,posit_y,posit_x,k);
  }
  else if( !w.is_empty() && g.is_empty())
  {
    if(s=="lda"){pred = miceFast::impute_raw(s,posit_y,posit_x,k);} else {pred = miceFast::imputeW(s,posit_y,posit_x,k);}
  }
  else if(!w.is_empty() && !g.is_empty())
  {
    if(s=="lda"){pred = miceFast::imputeby(s,posit_y,posit_x,k);} else {pred = miceFast::imputebyW(s,posit_y,posit_x,k);}
  }
  return pred;

}


arma::colvec miceFast::impute_raw(std::string s, int posit_y,arma::uvec posit_x,int k){

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  index_full = miceFast::get_index_full(posit_y, posit_x);
  index_NA = miceFast::get_index_NA(posit_y, posit_x);

  if(!x.col(posit_y).has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}

  arma::mat X_full = x(index_full,posit_x);
  arma::mat X_NA = x(index_NA,posit_x);
  arma::colvec Y_full = x(index_full,posit_y_uvec);

  arma::colvec pred = x(index_NA,posit_y_uvec);

  if(!(index_NA.n_elem==0) && ((index_full.n_elem>15 && s=="lda")|| (index_full.n_elem>posit_x.n_elem && s!="lda"))){

  pfunc f = funMap[s];
  pred = (*f)(Y_full,X_full,X_NA,k, ridge);

  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;

  return Y;
}


//Impute with grouping

arma::colvec miceFast::imputeby(std::string s, int posit_y,arma::uvec posit_x, int k){

  if(sorted == false){
    sortData_byg();
  }

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  index_full = miceFast::get_index_full(posit_y, posit_x);
  index_NA = miceFast::get_index_NA(posit_y, posit_x);

  //if(!Y.has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}

  if(!x.col(posit_y).has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(g.has_nan()){Rcpp::stop("There is NA values for the grouping variable");}

  arma::uvec g_int(N_rows);

  g_int = arma::conv_to<arma::uvec>::from(g);

  //grouping variable

  arma::uvec un = arma::unique(g_int);

  unsigned int group = un.n_elem;

  //quantitative model

  pfunc fun = funMap[s];

  //dividing data to NA and full

  arma::uvec g_full = g_int.elem(index_full);
  arma::uvec g_NA = g_int.elem(index_NA);

  //predictions container

  arma::colvec pred_all =  x(index_NA,posit_y_uvec) ;

  //iter

  arma::uvec his_full = arma::hist(g_full,un);
  arma::uvec ends_full = arma::cumsum(his_full);
  arma::uvec starts_full = arma::shift(ends_full,1) + 1;
  starts_full(0) = 1;

  arma::uvec his_NA = arma::hist(g_NA,un);
  arma::uvec ends_NA = arma::cumsum(his_NA);
  arma::uvec starts_NA = arma::shift(ends_NA,1) + 1;
  starts_NA(0) = 1;

  //unsigned int start = 0;
  //unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    int ss_NA = (int)starts_NA(a) - 1L;
    int ee_NA = (int)ends_NA(a) - 1L;
    int ss_full = (int)starts_full(a) - 1L;
    int ee_full = (int)ends_full(a) - 1L;

    if((ss_NA <= ee_NA) && (ss_full <= ee_full)){

      arma::mat X_full_0 = x(index_full.subvec(ss_full,ee_full),posit_x) ;
      arma::mat X_NA_0 = x(index_NA.subvec(ss_NA,ee_NA),posit_x);
      arma::colvec Y_full_0 =  x(index_full.subvec(ss_full,ee_full),posit_y_uvec) ;

      unsigned int N_obs_chunk = X_full_0.n_rows;

      if((N_obs_chunk<=15 && s=="lda")|| (N_obs_chunk<posit_x.n_elem && s!="lda")){ continue ;}

      arma::colvec pred =  (*fun)(Y_full_0,X_full_0,X_NA_0,k, ridge);

      //end = start + pred.n_elem - 1;

      pred_all.rows(ss_NA,ee_NA) = pred;

      //start = end + 1;

    }
  }

  arma::colvec Y = x.col(posit_y);

  Y.rows(index_NA) = pred_all;

  return Y;
}

//WEIGHTED

  typedef arma::colvec (*pfuncw)(arma::colvec&,arma::mat&,arma::colvec&,arma::mat&,int,double);
  std::map<std::string, pfuncw> funMapw = {{"lm_pred",fastLm_weighted},
  {"lm_noise",fastLm_weighted_noise},
  {"lm_bayes",fastLm_weighted_bayes},
  {"pmm",pmm_weighted_neibo}};


arma::colvec miceFast::imputeW(std::string s,int posit_y,arma::uvec posit_x,int k){

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  index_full = miceFast::get_index_full(posit_y, posit_x);
  index_NA = miceFast::get_index_NA(posit_y, posit_x);

  if(!x.col(posit_y).has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights");}
  if(arma::any(w<0)){Rcpp::stop("There are ngative values for the weights variable");}

  //dividing data to NA and full

  arma::mat X_full = x(index_full,posit_x);
  arma::mat X_NA = x(index_NA,posit_x);
  arma::colvec Y_full = x(index_full,posit_y_uvec);
  arma::colvec w_full = w.elem(index_full);

  arma::colvec pred = x(index_NA,posit_y_uvec);

  if(!(index_NA.n_elem==0) && ((X_full.n_rows>15 && s=="lda")|| (index_full.n_elem>posit_x.n_elem && s!="lda"))){

  pfuncw f = funMapw[s];
  pred = (*f)(Y_full,X_full,w_full,X_NA,k, ridge);

  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;

  return Y;
}

//Impute with grouping

arma::colvec miceFast::imputebyW(std::string s,int posit_y,arma::uvec posit_x,int k){

  if(!sorted){
    sortData_byg();
  }

  arma::uvec posit_y_uvec(1);
  posit_y_uvec(0) = posit_y;

  index_full = miceFast::get_index_full(posit_y, posit_x);
  index_NA = miceFast::get_index_NA(posit_y, posit_x);

  if(!x.col(posit_y).has_nan()){Rcpp::stop("There is no NA values for the dependent variable");}
  if(w.has_nan()){Rcpp::stop("There is NA values for weights variable");}
  if(g.has_nan()){Rcpp::stop("There is NA values for the grouping variable");}

  arma::uvec g_int(N_rows);

  g_int = arma::conv_to<arma::uvec>::from(g);

  //grouping variable

  arma::uvec un = arma::unique(g_int);
  unsigned int group = un.n_elem;

  //quantitative model

  pfuncw fun = funMapw[s];

  //dividing data to NA and full

  arma::uvec g_full = g_int.elem(index_full);
  arma::uvec g_NA = g_int.elem(index_NA);

  //predictions container

  arma::colvec pred_all =  x(index_NA,posit_y_uvec) ;

  // start end

  //iter

  arma::uvec his_full = arma::hist(g_full,un);
  arma::uvec ends_full = arma::cumsum(his_full);
  arma::uvec starts_full = arma::shift(ends_full,1) + 1;
  starts_full(0) = 1;

  arma::uvec his_NA = arma::hist(g_NA,un);
  arma::uvec ends_NA = arma::cumsum(his_NA);
  arma::uvec starts_NA = arma::shift(ends_NA,1) + 1;
  starts_NA(0) = 1;


  //unsigned int start = 0;
  //unsigned int end = 0;

  for(unsigned int  a=0;a<group;a++){

    int ss_NA = (int)starts_NA(a) - 1L;
    int ee_NA = (int)ends_NA(a) - 1L;
    int ss_full = (int)starts_full(a) - 1L;
    int ee_full = (int)ends_full(a) - 1L;

    if((ss_NA <= ee_NA) && (ss_full <= ee_full)){

    arma::mat X_full_0 = x(index_full.subvec(ss_full,ee_full),posit_x) ;
    arma::mat X_NA_0 = x(index_NA.subvec(ss_NA,ee_NA),posit_x);
    arma::colvec Y_full_0 =  x(index_full.subvec(ss_full,ee_full),posit_y_uvec) ;
    arma::colvec w_full_0 =  w(index_full.subvec(ss_full,ee_full));

    unsigned int N_obs_chunk = X_full_0.n_rows;

    if((N_obs_chunk<=15 && s=="lda")|| (N_obs_chunk<posit_x.n_elem && s!="lda")){ continue ;}

    arma::colvec pred =  (*fun)(Y_full_0,X_full_0,w_full_0,X_NA_0,k, ridge);

    //end = start + pred.n_elem - 1;

    pred_all.rows(ss_NA,ee_NA) = pred;

    //start = end + 1;

    }

  }

  arma::colvec Y = x.col(posit_y);

  Y.rows(index_NA) = pred_all;

  return Y;
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
    .method("get_ridge", &miceFast::get_ridge)
    .method("get_index", &miceFast::get_index)
    .method("is_sorted_byg", &miceFast::is_sorted_byg)
    .method("which_updated", &miceFast::which_updated)
    .method("vifs", &miceFast::vifs)
    .method("impute", &miceFast::impute)
    .method("impute_N", &miceFast::impute_N)
    .method("update_var", &miceFast::update_var)
    .method("get_models", &miceFast::get_models)
    .method("get_model", &miceFast::get_model)
    .method("set_data", &miceFast::set_data)
    .method("set_g", &miceFast::set_g)
    .method("set_w", &miceFast::set_w)
    .method("set_ridge", &miceFast::set_ridge)
    .method("sort_byg",&miceFast::sortData_byg)

  ;}
