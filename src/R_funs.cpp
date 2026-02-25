// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "miceFast.h"

// ---------------------------------------------------------------------------
// Index helpers (free-function versions for the R interface)
// ---------------------------------------------------------------------------

arma::uvec get_index_full_R(arma::mat &x, int posit_y, arma::uvec posit_x) {
  arma::colvec y_col = x.col(posit_y);
  arma::mat x_cols   = x.cols(posit_x);
  arma::uvec full_y  = complete_cases_vec(y_col);
  arma::uvec full_x  = complete_cases_mat(x_cols);
  return arma::find((full_y + full_x) == 2);
}

arma::uvec get_index_NA_R(arma::mat &x, int posit_y, arma::uvec posit_x) {
  arma::colvec y_col = x.col(posit_y);
  arma::mat x_cols   = x.cols(posit_x);
  arma::uvec full_y  = complete_cases_vec(y_col);
  arma::uvec full_x  = complete_cases_mat(x_cols);
  // y is missing AND x is complete
  return arma::find((1 - full_y) % full_x == 1);
}

// ---------------------------------------------------------------------------
// Minimum sample-size check shared by all imputation dispatchers
// ---------------------------------------------------------------------------
static bool has_enough_obs(arma::uword n_full, arma::uword n_x, const std::string &s) {
  if (s == "lda") return n_full > 15;
  return n_full > n_x;
}

// ---------------------------------------------------------------------------
// Core imputation dispatchers
// ---------------------------------------------------------------------------

// Unweighted dispatch
static arma::colvec impute_dispatch(arma::mat &x, const std::string &s,
                                     int posit_y, arma::uvec posit_x,
                                     int k, double ridge) {
  typedef arma::colvec (*pfunc)(arma::colvec &, arma::mat &, arma::mat &, int, double);
  static const std::map<std::string, pfunc> funMap = {
    {"lda",      fastLda},
    {"lm_pred",  fastLm_pred},
    {"lm_noise", fastLm_noise},
    {"lm_bayes", fastLm_bayes},
    {"pmm",      pmm_neibo}
  };

  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  arma::uvec idx_full = get_index_full_R(x, posit_y, posit_x);
  arma::uvec idx_NA   = get_index_NA_R(x, posit_y, posit_x);

  arma::colvec pred = x(idx_NA, posit_y_uvec);

  if (idx_NA.n_elem > 0 && has_enough_obs(idx_full.n_elem, posit_x.n_elem, s)) {
    arma::mat     X_full = x(idx_full, posit_x);
    arma::mat     X_NA   = x(idx_NA, posit_x);
    arma::colvec  Y_full = x(idx_full, posit_y_uvec);
    pfunc f = funMap.at(s);
    pred = f(Y_full, X_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(idx_NA) = pred;
  return Y;
}

// Weighted dispatch
static arma::colvec impute_dispatch_w(arma::mat &x, const std::string &s,
                                       int posit_y, arma::uvec posit_x,
                                       arma::colvec &w, int k, double ridge) {
  typedef arma::colvec (*pfuncw)(arma::colvec &, arma::mat &, arma::colvec &, arma::mat &, int, double);
  static const std::map<std::string, pfuncw> funMapw = {
    {"lm_pred",  fastLm_weighted},
    {"lm_noise", fastLm_weighted_noise},
    {"lm_bayes", fastLm_weighted_bayes},
    {"pmm",      pmm_weighted_neibo}
  };

  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  arma::uvec idx_full = get_index_full_R(x, posit_y, posit_x);
  arma::uvec idx_NA   = get_index_NA_R(x, posit_y, posit_x);

  arma::colvec pred = x(idx_NA, posit_y_uvec);

  if (idx_NA.n_elem > 0 && has_enough_obs(idx_full.n_elem, posit_x.n_elem, s)) {
    arma::mat     X_full = x(idx_full, posit_x);
    arma::mat     X_NA   = x(idx_NA, posit_x);
    arma::colvec  Y_full = x(idx_full, posit_y_uvec);
    arma::colvec  w_full = w.elem(idx_full);
    pfuncw f = funMapw.at(s);
    pred = f(Y_full, X_full, w_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(idx_NA) = pred;
  return Y;
}

// Kept for backward compatibility; delegates to the unified dispatchers above
arma::colvec impute_raw_R(arma::mat &x, std::string s, int posit_y,
                           arma::uvec posit_x, int k, double ridge) {
  return impute_dispatch(x, s, posit_y, posit_x, k, ridge);
}

arma::colvec imputeW_R(arma::mat &x, std::string s, int posit_y,
                        arma::uvec posit_x, arma::colvec w, int k, double ridge) {
  return impute_dispatch_w(x, s, posit_y, posit_x, w, k, ridge);
}

// ---------------------------------------------------------------------------
// VIF (Variance Inflation Factors)
// ---------------------------------------------------------------------------

static arma::mat cov2cor(const arma::mat &X) {
  arma::vec d = 1.0 / arma::sqrt(X.diag());
  return arma::diagmat(d) * X * arma::diagmat(d);
}

// [[Rcpp::export]]
arma::vec VIF_(arma::mat &x, int posit_y, arma::uvec posit_x,
               arma::uvec posit_x_var, bool correct) {

  arma::mat x_cols   = x.cols(posit_x - 1);
  arma::colvec y_col = x.col(posit_y - 1);
  arma::uvec Nvar_x_o = arma::unique(posit_x_var);
  int Nvar_o = Nvar_x_o.n_elem;

  arma::vec vifs(Nvar_o, arma::fill::none);

  arma::uvec full_rows = get_index_full_R(x, posit_y - 1, posit_x - 1);
  arma::mat full_x = x_cols.rows(full_rows);

  arma::uvec vols = arma::find(arma::var(full_x) > 0);
  int vols_n = vols.n_elem;

  if (vols_n < Nvar_o) {
    Rcpp::warning("There is at least a one zero variance variable");
    return vifs;
  }

  arma::mat x_vols = full_x.cols(vols);
  arma::uvec posit_x_var_v = posit_x_var.elem(vols);

  // Center columns
  x_vols.each_row() -= arma::mean(x_vols, 0);

  int Ncol = x_vols.n_cols;
  arma::uvec Nvar_x = arma::unique(posit_x_var_v);
  int Nvar = Nvar_x.n_elem;

  arma::mat xtx   = x_vols.t() * x_vols;
  arma::mat XXinv = arma::inv(xtx);
  if (Nvar < Ncol) {
    XXinv = cov2cor(XXinv);
  }
  double det_XXinv = arma::det(XXinv);

  arma::vec dfs(Nvar);
  arma::uvec pp = arma::linspace<arma::uvec>(0, Ncol - 1, Ncol);

  for (int i = 0; i < Nvar; i++) {
    int a = Nvar_x(i);
    arma::uvec ii = pp.elem(arma::find(posit_x_var_v == a));
    int ii_len  = ii.n_elem;
    int ii_first = ii(0);
    int ii_last  = ii(ii_len - 1);
    dfs(i) = ii_len;

    arma::mat XXinv_small = XXinv;
    XXinv_small.shed_rows(ii_first, ii_last);
    XXinv_small.shed_cols(ii_first, ii_last);

    vifs(i) = arma::as_scalar(arma::det(XXinv(ii, ii)) * arma::det(XXinv_small) / det_XXinv);
  }

  if (correct) {
    for (unsigned int i = 0; i < vifs.n_elem; i++) {
      vifs(i) = std::pow(vifs(i), 1.0 / (2.0 * dfs(i)));
    }
  }
  return vifs;
}

// ---------------------------------------------------------------------------
// R-facing fill_NA / fill_NA_N (unified logic)
// ---------------------------------------------------------------------------

// Common imputation logic for both fill_NA_ and fill_NA_N_
static arma::colvec fill_NA_impl(arma::mat &x, const std::string &model,
                                  int posit_y, arma::uvec posit_x,
                                  arma::colvec &w, int k, double ridge) {
  posit_x -= 1;
  posit_y -= 1;

  if (w.is_empty() || model == "lda") {
    return impute_dispatch(x, model, posit_y, posit_x, k, ridge);
  } else {
    return impute_dispatch_w(x, model, posit_y, posit_x, w, k, ridge);
  }
}

// [[Rcpp::export]]
arma::colvec fill_NA_N_(arma::mat &x, std::string model, int posit_y,
                         arma::uvec posit_x, arma::colvec w,
                         int k = 10, double ridge = 1e-6) {
  return fill_NA_impl(x, model, posit_y, posit_x, w, k, ridge);
}

// [[Rcpp::export]]
arma::colvec fill_NA_(arma::mat &x, std::string model, int posit_y,
                       arma::uvec posit_x, arma::colvec w,
                       double ridge = 1e-6) {
  return fill_NA_impl(x, model, posit_y, posit_x, w, 1, ridge);
}
