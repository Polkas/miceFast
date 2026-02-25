#ifndef MICEFAST_H
#define MICEFAST_H

#include <RcppArmadillo.h>
#include <string>
#include <map>

// =========================================================================
// miceFast class
// =========================================================================

class miceFast {
  arma::mat    x;          // data matrix
  arma::colvec g;          // grouping variable
  arma::colvec w;          // weight variable
  std::vector<int> updated;
  bool     sorted = false;
  unsigned int N_rows;
  unsigned int N_cols;
  double   ridge = 1e-6;

  arma::uvec index;
  arma::uvec index_NA;
  arma::uvec index_full;

public:
  miceFast();
  ~miceFast();

  // --- Imputation entry points ---
  Rcpp::List   impute(std::string s, int posit_y, arma::uvec posit_x);
  Rcpp::List   impute_N(std::string s, int posit_y, arma::uvec posit_x, int k);

  // --- Internal dispatch ---
  arma::colvec impute_raw(std::string s, int posit_y, arma::uvec posit_x, int k);
  arma::colvec imputeby(std::string s, int posit_y, arma::uvec posit_x, int k);
  arma::colvec imputeW(std::string s, int posit_y, arma::uvec posit_x, int k);
  arma::colvec imputebyW(std::string s, int posit_y, arma::uvec posit_x, int k);
  arma::colvec option_impute_multiple(std::string s, int posit_y, arma::uvec posit_x, int k);

  // --- Diagnostics ---
  arma::vec vifs(int posit_y, arma::uvec posit_x);
  std::string get_models(int posit_y);
  std::string get_model(int posit_y);

  // --- Index helpers ---
  arma::uvec get_index_full(int posit_y, arma::uvec posit_x);
  arma::uvec get_index_NA(int posit_y, arma::uvec posit_x);

  // --- Data management ---
  void sortData_byg();
  void set_data(arma::mat &_x);
  void set_g(arma::colvec &_g);
  void set_w(arma::colvec &_w);
  void set_ridge(double _ridge);
  void update_var(int posit_y, arma::vec impute);

  // --- Getters ---
  std::vector<int> which_updated();
  bool         is_sorted_byg();
  arma::mat    get_data();
  arma::colvec get_w();
  arma::colvec get_g();
  double       get_ridge();
  arma::uvec   get_index();
};

// =========================================================================
// Quantitative imputation models (unweighted)
// =========================================================================
arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);
arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);
arma::colvec fastLm_noise(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);
arma::colvec fastLda(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);
arma::colvec pmm_neibo(arma::colvec &y, arma::mat &X, arma::mat &X1, int k, double ridge);

// =========================================================================
// Quantitative imputation models (weighted)
// =========================================================================
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);
arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);
arma::colvec pmm_weighted_neibo(arma::colvec &y, arma::mat &X, arma::colvec &w, arma::mat &X1, int k, double ridge);

// =========================================================================
// Additional helper functions
// =========================================================================
arma::uvec complete_cases_mat(arma::mat &x);
arma::uvec complete_cases_vec(arma::colvec &y);
bool different_y_and_x(int posit_y, arma::uvec posit_x);
bool different_x(arma::uvec posit_x);
arma::mat sym(arma::mat x);

// =========================================================================
// R interface free functions
// =========================================================================
arma::uvec   get_index_full_R(arma::mat &x, int posit_y, arma::uvec posit_x);
arma::uvec   get_index_NA_R(arma::mat &x, int posit_y, arma::uvec posit_x);
arma::colvec impute_raw_R(arma::mat &x, std::string s, int posit_y, arma::uvec posit_x, int k, double ridge);
arma::colvec imputeW_R(arma::mat &x, std::string s, int posit_y, arma::uvec posit_x, arma::colvec w, int k, double ridge);
arma::vec    VIF(arma::mat &x, int posit_y, arma::uvec posit_x);
arma::colvec neibo(arma::colvec &y, arma::colvec &miss, int k);

#endif // MICEFAST_H
