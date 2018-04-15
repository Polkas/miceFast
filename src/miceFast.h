//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>

//Building class

class miceFast{

  arma::mat x;       //variables
  arma::colvec g;    //grouping
  arma::colvec w;//weights
  std::vector<int> updated;
  bool sorted = false;
  unsigned int N_rows;
  unsigned int N_cols;

  arma::uvec index;
  arma::uvec index_NA;
  arma::uvec index_full;

  public:

  miceFast();
  ~miceFast();

  Rcpp::List impute(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List impute_N(std::string s, int posit_y,arma::uvec posit_x,int times);
  arma::colvec impute_raw(std::string s, int posit_y,arma::uvec posit_x,int times);
  arma::colvec imputeby(std::string s, int posit_y,arma::uvec posit_x,int times);
  arma::colvec imputeW(std::string s, int posit_y,arma::uvec posit_x,int times);
  arma::colvec imputebyW(std::string s, int posit_y,arma::uvec posit_x,int times);
  arma::colvec option_impute_multiple(std::string s,int posit_y,arma::uvec posit_x,int times);
  arma::vec vifs(int posit_y,arma::uvec posit_x);

  std::string get_models(int posit_y);
  std::string get_model(int posit_y);
  arma::uvec get_index_full(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_NA(int posit_y, arma::uvec posit_x);

  void sortData_byg();
  void set_data(arma::mat& _x);
  void set_g(arma::colvec& _g);
  void set_w(arma::colvec& _w);
  void update_var(int posit_y,arma::vec impute);

  std::vector<int> which_updated();
  bool is_sorted_byg();
  arma::mat get_data();
  arma::colvec get_w();
  arma::colvec get_g();
  arma::uvec get_index();

};

//
//quantitative models
//

//Simple linear regression
arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1, int times);

//Weighted linear regression
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X,arma::colvec &w, arma::mat &X1, int times);

//weighted linear regression - noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1,int times);

//weighted linear regression - bayes
arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X,arma::colvec &w,arma::mat &X1, int times);

//linear regression - bayes
arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1, int times);

//Linear regression with noise
arma::colvec fastLm_noise(arma::colvec &y,arma::mat &X, arma::mat &X1, int times);

//LDA prediction model
arma::colvec fastLda( arma::colvec &y,  arma::mat &X, arma::mat &X1, int times);

//LDA prediction model - noise
//arma::colvec fastLda_noise( arma::colvec &y,  arma::mat &X, arma::mat &X1);

//QDA prediction model
//arma::colvec fastQda( arma::colvec &y,  arma::mat &X, arma::mat &X1);

//
//Additional functions
//

arma::uvec complete_cases_mat(arma::mat &x);

arma::uvec complete_cases_vec(arma::colvec &y);

arma::uvec histFast(arma::uvec &gg);

bool different_y_and_x(int posit_y, arma::uvec posit_x);

bool different_x(arma::uvec posit_x);
