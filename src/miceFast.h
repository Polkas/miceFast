
/*-------------------------------------------------------------------------------

This file is part of miceFast.

miceFast is free software: you can redistribute it and/or modify

it under the terms of the GNU General Public License as published by

the Free Software Foundation, either version 3 of the License, or

(at your option) any later version.


miceFast is distributed in the hope that it will be useful,

but WITHOUT ANY WARRANTY; without even the implied warranty of

MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the

GNU General Public License for more details.


You should have received a copy of the GNU General Public License

along with miceFast. If not, see <http://www.gnu.org/licenses/>.


Written by:

Maciej Nasinski

#-------------------------------------------------------------------------------*/

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
  arma::uvec index;

  public:

  miceFast();
  ~miceFast();

  Rcpp::List impute(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List impute_raw(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List imputeby(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List imputeW(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List imputebyW(std::string s, int posit_y,arma::uvec posit_x);
  Rcpp::List option_impute(std::string s,int posit_y,arma::uvec posit_x);
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
