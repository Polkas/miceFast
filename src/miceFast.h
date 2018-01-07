//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>

//Building class

class miceFast{

  arma::mat x;       //variables
  arma::uvec g;    //grouping
  bool sorted = false;    //sorted by g or not
  arma::colvec w;//weights
  bool updated = false;

public:

  miceFast(arma::mat _x);
  miceFast(arma::mat _x,arma::uvec _g,bool _s);
  miceFast(arma::mat _x,arma::colvec _w);
  miceFast(arma::mat _x,arma::uvec _g,bool _s,arma::colvec _w);

  Rcpp::List impute(std::string s, int posit_y,arma::uvec posit_x, bool force);
  arma::colvec impute_raw(std::string s, int posit_y,arma::uvec posit_x);
  arma::colvec imputeby(std::string s, int posit_y,arma::uvec posit_x);
  arma::colvec imputeW(std::string s, int posit_y,arma::uvec posit_x);
  arma::colvec imputebyW(std::string s, int posit_y,arma::uvec posit_x);

  std::string get_models(int posit_y);
  std::string get_model(int posit_y);
  arma::uvec get_index_full(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_NA(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_full_R(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_NA_R(int posit_y, arma::uvec posit_x);

  bool is_vars_updated();
  void sortData_byg();

};

//
//quantitative models
//

//Simple linear regression
arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1);

//Weighted linear regression
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X,arma::colvec &w, arma::mat &X1);

//weighted linear regression - noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1);

//weighted linear regression - bayes
arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1);

//linear regression - bayes
arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1);

//Linear regression with noise
arma::colvec fastLm_noise(arma::colvec &y,arma::mat &X, arma::mat &X1);

//LDA prediction model
arma::colvec fastLda( arma::colvec &y,  arma::mat &X, arma::mat &X1);

//
//Additional functions
//

arma::uvec complete_cases_mat(arma::mat &x);

arma::uvec complete_cases_vec(arma::colvec &y);

arma::uvec histFast(arma::uvec &gg);

bool different_y_and_x(int posit_y, arma::uvec posit_x);

bool different_x(arma::uvec posit_x);
