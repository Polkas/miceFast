#include <RcppArmadillo.h>
#include "miceFast.h"
#include <algorithm>
#include <vector>

// ---------------------------------------------------------------------------
// Internal helpers – DRY building blocks for the imputation models
// ---------------------------------------------------------------------------

// OLS fit result: coefficients, (X'X+ridge)^{-1}, residual sum-of-squares
struct OlsFit {
  arma::colvec coef;
  arma::mat    xxinv;      // (X'X + ridge*I)^{-1}
  double       rss;        // residual sum of squares
  int          n;          // number of observations
  int          p;          // number of predictors
};

// Fit unweighted OLS: coef = solve(X'X + ridge*I, X'y)
static OlsFit ols_fit(const arma::colvec &y, const arma::mat &X, double ridge) {
  OlsFit fit;
  fit.n = X.n_rows;
  fit.p = X.n_cols;

  arma::mat XX = X.t() * X;
  XX.diag() += ridge;

  fit.xxinv = arma::inv(XX);
  fit.coef  = fit.xxinv * (X.t() * y);

  arma::colvec res = y - X * fit.coef;
  fit.rss = arma::as_scalar(res.t() * res);
  return fit;
}

// Fit weighted OLS: transform y,X by sqrt(w), then standard OLS
static OlsFit wls_fit(const arma::colvec &y, const arma::mat &X,
                       const arma::colvec &w, double ridge) {
  arma::colvec wq = arma::sqrt(w);
  arma::colvec y2 = wq % y;
  arma::mat    X2 = X.each_col() % wq;

  // OLS on the sqrt(w)-transformed system gives:
  //   coef  = (X'WX + ridge*I)^{-1} X'Wy          (WLS solution)
  //   xxinv = (X'WX + ridge*I)^{-1}
  //   rss   = (y2-X2*coef)'(y2-X2*coef) = sum w_i*(y_i-x_i'coef)^2   (weighted RSS)
  // So ols_fit already produces the correct weighted RSS — no extra step needed.
  return ols_fit(y2, X2, ridge);
}

// Draw one set of Bayesian posterior coefficients + sigma
// Returns (perturbed_coef, sigma_draw)
struct BayesDraw {
  arma::colvec coef;
  double       sigma;
};

static BayesDraw posterior_draw(const OlsFit &fit) {
  BayesDraw draw;
  double df   = fit.n - fit.p;
  double chi2 = Rcpp::as<double>(Rcpp::rchisq(1, df));
  draw.sigma  = std::sqrt(fit.rss / chi2);

  arma::mat chol_xxinv = arma::chol(fit.xxinv);
  arma::colvec z(fit.p);
  z.randn();

  draw.coef = fit.coef + chol_xxinv.t() * z * draw.sigma;
  draw.coef.replace(arma::datum::nan, 0);
  return draw;
}

// Accumulate k stochastic predictions (noise model): X1 * coef + N(0, sigma)
static arma::colvec predict_noise_avg(const arma::colvec &coef, double sigma,
                                       const arma::mat &X1, int k) {
  int N_NA = X1.n_rows;
  arma::colvec pred_sum(N_NA, arma::fill::zeros);
  for (int i = 0; i < k; i++) {
    arma::vec noise(N_NA);
    noise.randn();
    pred_sum += X1 * coef + noise * sigma;
  }
  return pred_sum / static_cast<double>(k);
}

// Accumulate k Bayesian posterior-predictive draws
static arma::colvec predict_bayes_avg(const OlsFit &fit, const arma::mat &X1, int k) {
  int N_NA = X1.n_rows;
  arma::colvec pred_sum(N_NA, arma::fill::zeros);
  for (int i = 0; i < k; i++) {
    BayesDraw d = posterior_draw(fit);
    arma::vec noise(N_NA);
    noise.randn();
    pred_sum += X1 * d.coef + noise * d.sigma;
  }
  return pred_sum / static_cast<double>(k);
}

// ---------------------------------------------------------------------------
// Public model functions (unweighted)
// ---------------------------------------------------------------------------

// Simple deterministic linear regression prediction
arma::colvec fastLm_pred(arma::colvec &y, arma::mat &X, arma::mat &X1,
                          int /*k*/, double ridge) {
  OlsFit fit = ols_fit(y, X, ridge);
  return X1 * fit.coef;
}

// Linear regression with additive residual noise
arma::colvec fastLm_noise(arma::colvec &y, arma::mat &X, arma::mat &X1,
                           int k, double ridge) {
  OlsFit fit = ols_fit(y, X, ridge);
  double sigma = std::sqrt(fit.rss / (fit.n - fit.p));
  return predict_noise_avg(fit.coef, sigma, X1, k);
}

// Bayesian linear regression (posterior-predictive draws)
arma::colvec fastLm_bayes(arma::colvec &y, arma::mat &X, arma::mat &X1,
                           int k, double ridge) {
  OlsFit fit = ols_fit(y, X, ridge);
  return predict_bayes_avg(fit, X1, k);
}

// ---------------------------------------------------------------------------
// Public model functions (weighted)
// ---------------------------------------------------------------------------

// Weighted deterministic linear regression prediction
arma::colvec fastLm_weighted(arma::colvec &y, arma::mat &X, arma::colvec &w,
                              arma::mat &X1, int /*k*/, double ridge) {
  OlsFit fit = wls_fit(y, X, w, ridge);
  return X1 * fit.coef;
}

// Weighted linear regression with additive noise
arma::colvec fastLm_weighted_noise(arma::colvec &y, arma::mat &X, arma::colvec &w,
                                    arma::mat &X1, int k, double ridge) {
  OlsFit fit = wls_fit(y, X, w, ridge);
  double sigma = std::sqrt(fit.rss / (fit.n - fit.p));
  return predict_noise_avg(fit.coef, sigma, X1, k);
}

// Weighted Bayesian linear regression
arma::colvec fastLm_weighted_bayes(arma::colvec &y, arma::mat &X, arma::colvec &w,
                                    arma::mat &X1, int k, double ridge) {
  OlsFit fit = wls_fit(y, X, w, ridge);
  return predict_bayes_avg(fit, X1, k);
}

// ---------------------------------------------------------------------------
// LDA prediction model
// ---------------------------------------------------------------------------

arma::colvec fastLda(arma::colvec &y, arma::mat &X, arma::mat &X1,
                     int /*k*/, double ridge) {
  arma::uvec vars = arma::find(arma::var(X) > 0);
  arma::mat X_vol = X.cols(vars);

  const double tol = 1e-6;
  int N = X_vol.n_rows;
  int C = X_vol.n_cols;

  arma::vec un = arma::unique(y);
  int group = un.n_elem;
  if (group < 2 || group > 15) {
    Rcpp::stop("minimum 2 and maximum 15 categories");
  }

  arma::vec counts = arma::conv_to<arma::vec>::from(arma::hist(y, un));
  arma::vec prior  = counts / static_cast<double>(N);

  // Group means
  arma::mat group_means(group, C);
  for (int i = 0; i < group; i++) {
    arma::uvec idx = arma::find(y == un(i));
    group_means.row(i) = arma::mean(X_vol.rows(idx), 0);
  }

  // Expand group means to observation level
  arma::mat group_means_mat(N, C);
  for (int i = 0; i < N; i++) {
    group_means_mat.row(i) = group_means.rows(arma::find(un == y(i)));
  }

  // Within-group scaling
  arma::mat Sw = arma::sqrt(arma::cov(X_vol - group_means_mat));
  arma::mat scaling = arma::diagmat(1.0 / Sw.diag());

  double fac = 1.0 / static_cast<double>(N);
  arma::mat X0 = std::sqrt(fac) * (X_vol - group_means_mat) * scaling;

  arma::mat input = X0.t() * X0;
  input.diag() += ridge;

  arma::mat U;
  arma::vec s2;
  arma::mat V;
  arma::svd(U, s2, V, input);

  arma::vec s = arma::sqrt(s2);
  arma::uvec proper = arma::find(s > tol);
  int rank = proper.n_elem;

  if (rank == 0) {
    Rcpp::stop("rank = 0: variables are numerically constant");
  }
  if (rank < C) {
    Rcpp::warning("variables are collinear");
  }

  scaling = scaling * U.cols(proper) * arma::diagmat(1.0 / s.elem(proper));

  arma::vec xb = group_means.t() * prior;
  fac = 1.0 / static_cast<double>(C);

  arma::mat X_s = arma::conv_to<arma::rowvec>::from(
    arma::sqrt(static_cast<double>(N) * prior * fac)) *
    (group_means.each_row() - xb.t()) * scaling;

  arma::svd(U, s2, V, X_s.t() * X_s);
  s = arma::sqrt(s2);
  proper = arma::find(s > tol * s(0));
  rank = proper.n_elem;

  if (rank == 0) {
    Rcpp::stop("group means are numerically identical");
  }

  scaling = scaling * U.cols(proper);

  arma::vec dist_base = 0.5 * arma::mean(arma::pow(group_means * scaling, 2), 1)
                        - arma::log(prior);

  // Prediction
  int N_pred = X1.n_rows;
  arma::mat X1_vol = X1.cols(vars);
  arma::mat dist_raw = arma::repmat(dist_base.t(), N_pred, 1)
                      - X1_vol * scaling * (group_means * scaling).t();
  dist_raw.each_col() -= arma::min(dist_raw, 1);
  dist_raw.transform([](double val) { return -std::exp(val); });

  arma::colvec pred(N_pred);
  for (int i = 0; i < N_pred; i++) {
    pred(i) = un(arma::index_max(dist_raw.row(i)));
  }
  return pred;
}

// ---------------------------------------------------------------------------
// Predictive Mean Matching (PMM)
// ---------------------------------------------------------------------------

//' Finding in random manner one of the k closest points in a certain vector
//' for each value in a second vector
//'
//' @description This function uses pre-sorting of y and binary search to find
//'   one of the k closest values for each miss.
//'
//' @param y numeric vector values to be looked up
//' @param miss numeric vector values to be looked for
//' @param k integer number of nearest neighbours to sample from
//'
//' @return a numeric vector
//'
//' @name neibo
//' @export
// [[Rcpp::export]]
arma::colvec neibo(arma::colvec &y, arma::colvec &miss, int k) {
  int n_y = y.n_rows;
  int n_miss = miss.n_rows;

  k = std::min(k, n_y);
  k = std::max(k, 1);

  arma::colvec result(n_miss);
  arma::vec which_n = arma::floor(Rcpp::as<arma::vec>(Rcpp::runif(n_miss, 0, k)));
  arma::uvec which = arma::conv_to<arma::uvec>::from(which_n);

  std::vector<double> y_sorted = arma::conv_to<std::vector<double>>::from(y);
  std::sort(y_sorted.begin(), y_sorted.end());

  for (int i = 0; i < n_miss; i++) {
    double mm = miss[i];
    int count = 0;
    std::vector<int> resus(k);

    auto iter_geq = std::lower_bound(y_sorted.begin(), y_sorted.end(), mm);
    int r = iter_geq - y_sorted.begin();
    int l = r - 1;

    // Expand outward from the insertion point, picking the closer side
    while (l >= 0 && r < n_y && count < k) {
      if (mm - y_sorted[l] < y_sorted[r] - mm)
        resus[count++] = l--;
      else
        resus[count++] = r++;
    }
    while (count < k && l >= 0)
      resus[count++] = l--;
    while (count < k && r < n_y)
      resus[count++] = r++;

    result[i] = y_sorted[resus[which[i]]];
  }
  return result;
}

// PMM core: fit Bayesian model, draw predictions, match to observed
static arma::colvec pmm_core(const OlsFit &fit, const arma::colvec &y,
                              const arma::mat &X, const arma::mat &X1,
                              int k) {
  BayesDraw d = posterior_draw(fit);

  int N_NA = X1.n_rows;
  arma::vec noise(N_NA);
  noise.randn();

  arma::colvec ypred_mis  = X1 * d.coef + noise * d.sigma;
  arma::colvec ypred_full = X * fit.coef;
  return neibo(ypred_full, ypred_mis, k);
}

// [[Rcpp::export]]
arma::colvec pmm_weighted_neibo(arma::colvec &y, arma::mat &X, arma::colvec &w,
                                 arma::mat &X1, int k, double ridge) {
  OlsFit fit = wls_fit(y, X, w, ridge);
  return pmm_core(fit, y, X, X1, k);
}

// [[Rcpp::export]]
arma::colvec pmm_neibo(arma::colvec &y, arma::mat &X, arma::mat &X1,
                        int k, double ridge) {
  OlsFit fit = ols_fit(y, X, ridge);
  return pmm_core(fit, y, X, X1, k);
}
