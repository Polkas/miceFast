#include <RcppArmadillo.h>
#include <random>

class corrData {
  int nr_cat;
  int n_row;
  int n_col;
  arma::vec mu;
  arma::mat cors;

public:
  corrData(int _n_row, arma::vec _mu, arma::mat _cors)
    : nr_cat(2), n_row(_n_row), n_col(_mu.n_elem), mu(_mu), cors(_cors) {}

  corrData(int _nr_cat, int _n_row, arma::vec _mu, arma::mat _cors)
    : nr_cat(_nr_cat), n_row(_n_row), n_col(_mu.n_elem), mu(_mu), cors(_cors) {}

  arma::mat fill(std::string type);
};

arma::mat corrData::fill(std::string type) {
  // Generate standard normal random data
  arma::mat X(n_row, n_col);
  std::mt19937 engine;
  std::normal_distribution<double> distr(0.0, 1.0);
  X.imbue([&]() { return distr(engine); });

  // Apply correlation structure
  arma::mat X_new = X * arma::inv(arma::chol(arma::cor(X))) * arma::chol(cors);

  // Discretize first column if needed
  int n_categories = 0;
  if (type == "binom") {
    n_categories = 2;
  } else if (type == "discrete") {
    n_categories = nr_cat;
  } else if (type != "contin") {
    Rcpp::stop("type must be one of: binom, discrete, contin");
  }

  if (n_categories > 0) {
    arma::colvec col0 = X_new.col(0);
    for (int i = 0; i < n_row; i++) {
      col0(i) = Rcpp::stats::pnorm_0(col0(i), 1, 0) * n_categories;
    }
    X_new.col(0) = arma::ceil(col0);
  }

  // Add means to all columns except the first (which was discretized)
  for (int i = 1; i < n_col; i++) {
    X_new.col(i) += mu(i);
  }

  return X_new;
}

RCPP_MODULE(corrData) {
  using namespace Rcpp;
  class_<corrData>("corrData")
    .constructor<int, arma::vec, arma::mat>()
    .constructor<int, int, arma::vec, arma::mat>()
    .method("fill", &corrData::fill);
}
