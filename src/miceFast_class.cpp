#include <RcppArmadillo.h>
#include <string>
#include <map>
#include "miceFast.h"

// =========================================================================
// Constructor / Destructor
// =========================================================================

miceFast::miceFast() {}
miceFast::~miceFast() {}

// =========================================================================
// Data management
// =========================================================================

void miceFast::set_data(arma::mat &_x) {
  x = arma::mat(_x.begin(), _x.n_rows, _x.n_cols, false, false);
  N_rows = x.n_rows;
  N_cols = x.n_cols;
  index = arma::regspace<arma::uvec>(1, N_rows);
}

void miceFast::set_g(arma::colvec &_g) {
  if (x.is_empty()) Rcpp::stop("There is no data provided");
  if (N_rows != _g.n_rows) Rcpp::stop("Wrong number of elements");
  if (_g.has_nan()) Rcpp::stop("There are NA values for the grouping variable");
  g = arma::colvec(_g.begin(), _g.n_rows, false, false);
  sorted = g.is_sorted();
}

void miceFast::set_w(arma::colvec &_w) {
  if (x.is_empty()) Rcpp::stop("There is no data provided");
  if (N_rows != _w.n_rows) Rcpp::stop("Wrong number of elements");
  if (_w.has_nan()) Rcpp::stop("There are NA values for weights");
  if (arma::any(_w < 0)) Rcpp::stop("There are negative values for the weights variable");
  w = arma::colvec(_w.begin(), _w.n_rows, false, false);
}

void miceFast::set_ridge(double _ridge) { ridge = _ridge; }

void miceFast::update_var(int posit_y, arma::vec impute) {
  if (x.is_empty()) Rcpp::stop("at least set the data");
  if (N_rows != impute.n_elem) Rcpp::stop("wrong number of observations");
  int col = posit_y - 1;
  x.col(col) = impute;
  updated.push_back(posit_y);
}

void miceFast::sortData_byg() {
  if (g.is_empty()) Rcpp::stop("There is no a grouping variable provided");
  if (!sorted) {
    Rcpp::warning("\n Data was sorted by the grouping variable - use `get_index()` to retrieve an order");
    arma::uvec order = arma::stable_sort_index(g);
    x = x.rows(order);
    g = g.rows(order);
    index = index.elem(order);
    if (!w.is_empty()) {
      w = w.rows(order);
    }
    sorted = true;
  }
}

// =========================================================================
// Getters
// =========================================================================

bool miceFast::is_sorted_byg() {
  if (g.is_empty()) Rcpp::stop("There is no grouping variable provided");
  return sorted;
}

arma::mat    miceFast::get_data()  { if (x.is_empty()) Rcpp::stop("There is no data provided"); return x; }
arma::colvec miceFast::get_w()     { if (w.is_empty()) Rcpp::stop("There is no weighting variable provided"); return w; }
arma::colvec miceFast::get_g()     { if (g.is_empty()) Rcpp::stop("There is no grouping variable provided"); return g; }
double       miceFast::get_ridge() { return ridge; }
arma::uvec   miceFast::get_index() { if (x.is_empty()) Rcpp::stop("There is no data provided"); return index; }
std::vector<int> miceFast::which_updated() { return updated; }

// =========================================================================
// Index helpers
// =========================================================================

arma::uvec miceFast::get_index_full(int posit_y, arma::uvec posit_x) {
  arma::colvec y_col = x.col(posit_y);
  arma::mat    x_cols = x.cols(posit_x);
  arma::uvec full_y = complete_cases_vec(y_col);
  arma::uvec full_x = complete_cases_mat(x_cols);
  return arma::find((full_y + full_x) == 2);
}

arma::uvec miceFast::get_index_NA(int posit_y, arma::uvec posit_x) {
  arma::colvec y_col = x.col(posit_y);
  arma::mat    x_cols = x.cols(posit_x);
  arma::uvec full_y = complete_cases_vec(y_col);
  arma::uvec full_x = complete_cases_mat(x_cols);
  // y is missing AND x is complete
  return arma::find((1 - full_y) % full_x == 1);
}

// =========================================================================
// Diagnostics
// =========================================================================

arma::vec miceFast::vifs(int posit_y, arma::uvec posit_x) {
  if (!different_y_and_x(posit_y, posit_x)) Rcpp::stop("the same variable is dependent and independent");
  if (!different_x(posit_x)) Rcpp::stop("the same variables repeated few times as independent");
  if (x.is_empty()) Rcpp::stop("at least set the data");

  arma::mat x_cols = x.cols(posit_x - 1);
  arma::uvec full_rows = get_index_full(posit_y - 1, posit_x - 1);
  arma::mat full_x = x_cols.rows(full_rows);

  // Center columns
  full_x.each_row() -= arma::mean(full_x, 0);

  arma::mat XXinv = arma::inv(full_x.t() * full_x);
  int Ncol = full_x.n_cols;
  double det_XXinv = arma::det(XXinv);

  arma::vec vifs_out(Ncol);
  for (int i = 0; i < Ncol; i++) {
    arma::mat XXinv_small = XXinv;
    XXinv_small.shed_row(i);
    XXinv_small.shed_col(i);
    vifs_out(i) = XXinv(i, i) * arma::det(XXinv_small) / det_XXinv;
  }
  return vifs_out;
}

// Helper: classify a variable by its unique-value count
static std::string classify_variable(const arma::colvec &Y_raw, bool verbose) {
  arma::uvec idx_full = arma::find_finite(Y_raw);
  if (!Y_raw.has_nan()) return "no NA values for the dependent variable";

  arma::colvec Y = Y_raw.rows(idx_full);
  arma::colvec uY = arma::unique(Y);
  int un_n = uY.n_elem;

  if (un_n <= 1) return verbose ? "one unique value" : "one unique value";
  if (un_n == 2) return verbose ? "recommended lda or (lm_pred,lm_bayes,lm_noise, pmm - round results)" : "lda";
  if (un_n <= 15) return verbose ? "recommended lda or (lm_pred,lm_bayes,lm_noise, pmm - remember to round results if needed)" : "lda";
  return verbose ? "lm_pred or lm_bayes or lm_noise or pmm" : "lm_pred";
}

std::string miceFast::get_models(int posit_y) {
  if (x.is_empty()) Rcpp::stop("at least set the data");
  return classify_variable(x.col(posit_y - 1), true);
}

std::string miceFast::get_model(int posit_y) {
  if (x.is_empty()) Rcpp::stop("at least set the data");
  return classify_variable(x.col(posit_y - 1), false);
}

// =========================================================================
// Imputation – common validation + return logic
// =========================================================================

static void validate_impute_args(int posit_y, arma::uvec posit_x, const arma::mat &x) {
  if (!different_y_and_x(posit_y, posit_x)) Rcpp::stop("the same variable is dependent and independent");
  if (!different_x(posit_x)) Rcpp::stop("the same variables repeated few times as independent");
  if (x.is_empty()) Rcpp::stop("at least set the data");
}

static Rcpp::List make_impute_result(const arma::colvec &pred,
                                      const arma::uvec &idx_NA,
                                      const arma::uvec &idx_full,
                                      arma::uword N) {
  arma::uvec index_NA_ret(N, arma::fill::zeros);
  index_NA_ret.elem(idx_NA).fill(1);
  arma::uvec index_full_ret(N, arma::fill::zeros);
  index_full_ret.elem(idx_full).fill(1);
  return Rcpp::List::create(
    Rcpp::Named("imputations")   = pred,
    Rcpp::Named("index_imputed") = index_NA_ret,
    Rcpp::Named("index_full")    = index_full_ret
  );
}

Rcpp::List miceFast::impute(std::string s, int posit_y, arma::uvec posit_x) {
  validate_impute_args(posit_y, posit_x, x);
  posit_x -= 1;
  posit_y -= 1;
  arma::colvec pred = option_impute_multiple(s, posit_y, posit_x, 1);
  return make_impute_result(pred, index_NA, index_full, x.n_rows);
}

Rcpp::List miceFast::impute_N(std::string s, int posit_y, arma::uvec posit_x, int k) {
  if (s != "lm_bayes" && s != "lm_noise" && s != "pmm") {
    Rcpp::stop("Works only for `lm_bayes`, `lm_noise` and `pmm` models");
  }
  validate_impute_args(posit_y, posit_x, x);
  posit_x -= 1;
  posit_y -= 1;
  arma::colvec pred = option_impute_multiple(s, posit_y, posit_x, k);
  return make_impute_result(pred, index_NA, index_full, x.n_rows);
}

// =========================================================================
// Imputation dispatch
// =========================================================================

// Function pointer maps (unweighted / weighted)
typedef arma::colvec (*pfunc)(arma::colvec &, arma::mat &, arma::mat &, int, double);
static const std::map<std::string, pfunc> funMap = {
  {"lda",      fastLda},
  {"lm_pred",  fastLm_pred},
  {"lm_noise", fastLm_noise},
  {"lm_bayes", fastLm_bayes},
  {"pmm",      pmm_neibo}
};

typedef arma::colvec (*pfuncw)(arma::colvec &, arma::mat &, arma::colvec &, arma::mat &, int, double);
static const std::map<std::string, pfuncw> funMapw = {
  {"lm_pred",  fastLm_weighted},
  {"lm_noise", fastLm_weighted_noise},
  {"lm_bayes", fastLm_weighted_bayes},
  {"pmm",      pmm_weighted_neibo}
};

// Minimum sample-size guard (same logic as R_funs.cpp but inline here)
static bool enough_obs(arma::uword n_full, arma::uword n_x, const std::string &s) {
  return (s == "lda") ? n_full > 15 : n_full > n_x;
}

arma::colvec miceFast::option_impute_multiple(std::string s, int posit_y,
                                               arma::uvec posit_x, int k) {
  bool has_w = !w.is_empty();
  bool has_g = !g.is_empty();
  bool use_w = has_w && (s != "lda");  // LDA ignores weights

  if (has_g) {
    return use_w ? imputebyW(s, posit_y, posit_x, k)
                 : imputeby(s, posit_y, posit_x, k);
  } else {
    return use_w ? imputeW(s, posit_y, posit_x, k)
                 : impute_raw(s, posit_y, posit_x, k);
  }
}

// =========================================================================
// Non-grouped imputation
// =========================================================================

arma::colvec miceFast::impute_raw(std::string s, int posit_y,
                                   arma::uvec posit_x, int k) {
  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  index_full = get_index_full(posit_y, posit_x);
  index_NA   = get_index_NA(posit_y, posit_x);

  if (!x.col(posit_y).has_nan()) Rcpp::stop("There is no NA values for the dependent variable");

  arma::colvec pred = x(index_NA, posit_y_uvec);

  if (index_NA.n_elem > 0 && enough_obs(index_full.n_elem, posit_x.n_elem, s)) {
    arma::mat    X_full = x(index_full, posit_x);
    arma::mat    X_NA   = x(index_NA, posit_x);
    arma::colvec Y_full = x(index_full, posit_y_uvec);
    pred = funMap.at(s)(Y_full, X_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;
  return Y;
}

arma::colvec miceFast::imputeW(std::string s, int posit_y,
                                arma::uvec posit_x, int k) {
  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  index_full = get_index_full(posit_y, posit_x);
  index_NA   = get_index_NA(posit_y, posit_x);

  if (!x.col(posit_y).has_nan()) Rcpp::stop("There are no NA values for the dependent variable");

  arma::colvec pred = x(index_NA, posit_y_uvec);

  if (index_NA.n_elem > 0 && enough_obs(index_full.n_elem, posit_x.n_elem, s)) {
    arma::mat    X_full = x(index_full, posit_x);
    arma::mat    X_NA   = x(index_NA, posit_x);
    arma::colvec Y_full = x(index_full, posit_y_uvec);
    arma::colvec w_full = w.elem(index_full);
    pred = funMapw.at(s)(Y_full, X_full, w_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred;
  return Y;
}

// =========================================================================
// Grouped imputation – shared loop scaffolding
// =========================================================================

// Helper struct to hold group-range bookkeeping
struct GroupRanges {
  arma::uvec un;
  arma::uvec starts_full, ends_full;
  arma::uvec starts_NA,   ends_NA;
  unsigned int n_groups;
};

static GroupRanges compute_group_ranges(const arma::colvec &g,
                                         const arma::uvec &idx_full,
                                         const arma::uvec &idx_NA,
                                         unsigned int N_rows) {
  GroupRanges gr;
  arma::uvec g_int = arma::conv_to<arma::uvec>::from(g);
  gr.un = arma::unique(g_int);
  gr.n_groups = gr.un.n_elem;

  arma::uvec g_full = g_int.elem(idx_full);
  arma::uvec his_full = arma::hist(g_full, gr.un);
  gr.ends_full   = arma::cumsum(his_full);
  gr.starts_full = arma::shift(gr.ends_full, 1) + 1;
  gr.starts_full(0) = 1;

  arma::uvec g_NA = g_int.elem(idx_NA);
  arma::uvec his_NA = arma::hist(g_NA, gr.un);
  gr.ends_NA   = arma::cumsum(his_NA);
  gr.starts_NA = arma::shift(gr.ends_NA, 1) + 1;
  gr.starts_NA(0) = 1;

  return gr;
}

arma::colvec miceFast::imputeby(std::string s, int posit_y,
                                 arma::uvec posit_x, int k) {
  if (!sorted) sortData_byg();

  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  index_full = get_index_full(posit_y, posit_x);
  index_NA   = get_index_NA(posit_y, posit_x);

  if (!x.col(posit_y).has_nan()) Rcpp::stop("There are no NA values for the dependent variable");

  GroupRanges gr = compute_group_ranges(g, index_full, index_NA, N_rows);
  pfunc fun = funMap.at(s);
  arma::colvec pred_all = x(index_NA, posit_y_uvec);

  for (unsigned int a = 0; a < gr.n_groups; a++) {
    int ss_NA   = static_cast<int>(gr.starts_NA(a)) - 1;
    int ee_NA   = static_cast<int>(gr.ends_NA(a)) - 1;
    int ss_full = static_cast<int>(gr.starts_full(a)) - 1;
    int ee_full = static_cast<int>(gr.ends_full(a)) - 1;

    if (ss_NA > ee_NA || ss_full > ee_full) continue;

    arma::mat    X_full = x(index_full.subvec(ss_full, ee_full), posit_x);
    arma::mat    X_NA   = x(index_NA.subvec(ss_NA, ee_NA), posit_x);
    arma::colvec Y_full = x(index_full.subvec(ss_full, ee_full), posit_y_uvec);

    if (!enough_obs(X_full.n_rows, posit_x.n_elem, s)) continue;

    pred_all.rows(ss_NA, ee_NA) = fun(Y_full, X_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred_all;
  return Y;
}

arma::colvec miceFast::imputebyW(std::string s, int posit_y,
                                  arma::uvec posit_x, int k) {
  if (!sorted) sortData_byg();

  arma::uvec posit_y_uvec = {static_cast<arma::uword>(posit_y)};
  index_full = get_index_full(posit_y, posit_x);
  index_NA   = get_index_NA(posit_y, posit_x);

  if (!x.col(posit_y).has_nan()) Rcpp::stop("There are no NA values for the dependent variable");

  GroupRanges gr = compute_group_ranges(g, index_full, index_NA, N_rows);
  pfuncw fun = funMapw.at(s);
  arma::colvec pred_all = x(index_NA, posit_y_uvec);

  for (unsigned int a = 0; a < gr.n_groups; a++) {
    int ss_NA   = static_cast<int>(gr.starts_NA(a)) - 1;
    int ee_NA   = static_cast<int>(gr.ends_NA(a)) - 1;
    int ss_full = static_cast<int>(gr.starts_full(a)) - 1;
    int ee_full = static_cast<int>(gr.ends_full(a)) - 1;

    if (ss_NA > ee_NA || ss_full > ee_full) continue;

    arma::mat    X_full = x(index_full.subvec(ss_full, ee_full), posit_x);
    arma::mat    X_NA   = x(index_NA.subvec(ss_NA, ee_NA), posit_x);
    arma::colvec Y_full = x(index_full.subvec(ss_full, ee_full), posit_y_uvec);
    arma::colvec w_full = w(index_full.subvec(ss_full, ee_full));

    if (!enough_obs(X_full.n_rows, posit_x.n_elem, s)) continue;

    pred_all.rows(ss_NA, ee_NA) = fun(Y_full, X_full, w_full, X_NA, k, ridge);
  }

  arma::colvec Y = x.col(posit_y);
  Y.rows(index_NA) = pred_all;
  return Y;
}

// =========================================================================
// RCPP Module registration
// =========================================================================

RCPP_MODULE(miceFast) {
  using namespace Rcpp;
  class_<miceFast>("miceFast")
    .default_constructor()
    .method("get_data",       &miceFast::get_data)
    .method("get_w",          &miceFast::get_w)
    .method("get_g",          &miceFast::get_g)
    .method("get_ridge",      &miceFast::get_ridge)
    .method("get_index",      &miceFast::get_index)
    .method("is_sorted_byg",  &miceFast::is_sorted_byg)
    .method("which_updated",  &miceFast::which_updated)
    .method("vifs",           &miceFast::vifs)
    .method("impute",         &miceFast::impute)
    .method("impute_N",       &miceFast::impute_N)
    .method("update_var",     &miceFast::update_var)
    .method("get_models",     &miceFast::get_models)
    .method("get_model",      &miceFast::get_model)
    .method("set_data",       &miceFast::set_data)
    .method("set_g",          &miceFast::set_g)
    .method("set_w",          &miceFast::set_w)
    .method("set_ridge",      &miceFast::set_ridge)
    .method("sort_byg",       &miceFast::sortData_byg);
}
