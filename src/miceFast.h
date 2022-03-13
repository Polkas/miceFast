//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <string>

//Building class

class miceFast{

  arma::mat x;       //variables
  arma::sp_mat x_sp;
  arma::colvec g;    //grouping
  arma::colvec w;//weights
  std::vector<int> updated;
  bool sorted = false;
  unsigned int N_rows;
  unsigned int N_cols;
  double ridge = 1e-6;

  arma::uvec index;
  arma::uvec index_NA;
  arma::uvec index_full;

  public:

  miceFast();
  ~miceFast();

  Rcpp::List impute(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List impute_N(std::string s, int posit_y,arma::uvec posit_x,int k);
  arma::colvec impute_raw(std::string s, int posit_y,arma::uvec posit_x,int k);
  arma::colvec imputeby(std::string s, int posit_y,arma::uvec posit_x,int k);
  arma::colvec imputeW(std::string s, int posit_y,arma::uvec posit_x,int k);
  arma::colvec imputebyW(std::string s, int posit_y,arma::uvec posit_x,int k);
  arma::colvec option_impute_multiple(std::string s,int posit_y,arma::uvec posit_x,int k);
  arma::vec vifs(int posit_y,arma::uvec posit_x);

  std::string get_models(int posit_y);
  std::string get_model(int posit_y);
  arma::uvec get_index_full(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_NA(int posit_y, arma::uvec posit_x);

  void sortData_byg();
  void set_data(arma::mat& _x);
  void set_data_sparse(arma::sp_mat & _x);
  void set_g(arma::colvec& _g);
  void set_w(arma::colvec& _w);
  void set_ridge(double _ridge);
  void update_var(int posit_y,arma::vec impute);

  std::vector<int> which_updated();
  bool is_sorted_byg();
  arma::mat get_data();
  arma::colvec get_w();
  arma::colvec get_g();
  double get_ridge();
  arma::uvec get_index();

};


const double ridge = 0.00001;

//
//quantitative models
//

//Simple linear regression
arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);

//Weighted linear regression
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);

//weighted linear regression - noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1,int k, double ridge);

//weighted linear regression - bayes
arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);

//linear regression - bayes
arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);

//Linear regression with noise
arma::colvec fastLm_noise(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);

//LDA prediction model
arma::colvec fastLda(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);

//PMM
arma::colvec pmm_weighted_neibo(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);

arma::colvec pmm_neibo( arma::colvec &y, arma::mat &X,arma::mat &X1,int k, double ridge );

//LDA prediction model - noise
//arma::colvec fastLda_noise( arma::colvec &y,  arma::mat &X, arma::mat &X1);

//QDA prediction model
//arma::colvec fastQda( arma::colvec &y,  arma::mat &X, arma::mat &X1);

//
//Additional functions
//

arma::uvec complete_cases_mat(arma::mat &x);

arma::uvec complete_cases_vec(arma::colvec &y);

bool different_y_and_x(int posit_y, arma::uvec posit_x);

bool different_x(arma::uvec posit_x);

arma::mat sym(arma::mat x);
//
//R interface
//

arma::uvec get_index_full_R(arma::mat &x, int posit_y, arma::uvec posit_x);

arma::uvec get_index_NA_R(arma::mat &x, int posit_y, arma::uvec posit_x);

arma::colvec impute_raw_R(arma::mat &x, std::string s, int posit_y, arma::uvec posit_x, int k, double ridge);

arma::colvec imputeW_R(arma::mat &x, std::string s, int posit_y, arma::uvec posit_x, arma::colvec w,int k, double ridge);

arma::vec VIF(arma::mat &x,int posit_y, arma::uvec posit_x);

arma::colvec fill_NA_N(arma::mat &x, std::string model, int posit_y,arma::uvec posit_x, arma::colvec w, int k, double ridge);

arma::colvec fill_NA(arma::mat &x, std::string model, int posit_y, arma::uvec posit_x, arma::colvec w, double ridge);

arma::colvec neibo(arma::colvec &y, arma::colvec &miss, int k);
