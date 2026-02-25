#' Pool results from models fitted on multiply imputed datasets
#'
#' @description
#' Combines parameter estimates and standard errors from models fitted on
#' \emph{m} multiply imputed datasets using Rubin's rules (Rubin, 1987).
#' Degrees of freedom are adjusted using the Barnard-Rubin (1999) small-sample correction.
#'
#' This function works with any fitted model that supports \code{\link[stats]{coef}} and
#' \code{\link[stats]{vcov}} methods (e.g., \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' \code{survival::coxph}, etc.).
#'
#' Results are validated against \code{\link[mice]{pool}} from the \pkg{mice} package
#' for \code{lm}, \code{glm} (logistic and Poisson), weighted regression, interactions,
#' and varying numbers of imputations.
#'
#' @param fits a list of fitted model objects of length \emph{m} >= 2.
#'   Each model must support \code{coef()} and \code{vcov()} methods.
#'   All models must have the same number of coefficients.
#' @param dfcom a positive integer or \code{Inf}. The complete-data degrees of freedom.
#'   If \code{NULL} (default), it is extracted from the fitted models via
#'   \code{\link[stats]{df.residual}}. Set to \code{Inf} to skip the Barnard-Rubin
#'   small-sample correction.
#'
#' @return A data.frame with one row per parameter and columns:
#'   \describe{
#'     \item{term}{Coefficient name.}
#'     \item{m}{Number of imputations.}
#'     \item{estimate}{Pooled estimate (average across \emph{m} models).}
#'     \item{std.error}{Pooled standard error (\code{sqrt(t)}).}
#'     \item{statistic}{t-statistic (\code{estimate / std.error}).}
#'     \item{p.value}{Two-sided p-value from a t-distribution with \code{df} degrees of freedom.}
#'     \item{df}{Degrees of freedom (Barnard-Rubin adjusted).}
#'     \item{riv}{Relative increase in variance due to nonresponse: \code{(1 + 1/m) * b / ubar}.}
#'     \item{lambda}{Proportion of total variance attributable to missingness: \code{(1 + 1/m) * b / t}.}
#'     \item{fmi}{Fraction of missing information.}
#'     \item{ubar}{Within-imputation variance (average of the \emph{m} variance estimates).}
#'     \item{b}{Between-imputation variance (variance of the \emph{m} point estimates).}
#'     \item{t}{Total variance: \code{ubar + (1 + 1/m) * b}.}
#'     \item{dfcom}{Complete-data degrees of freedom used.}
#'     \item{conf.low}{Lower bound of the 95\% confidence interval.}
#'     \item{conf.high}{Upper bound of the 95\% confidence interval.}
#'   }
#'
#' @references
#' Rubin, D.B. (1987). \emph{Multiple Imputation for Nonresponse in Surveys}.
#' John Wiley & Sons.
#'
#' Barnard, J. and Rubin, D.B. (1999). Small-sample degrees of freedom with
#' multiple imputation. \emph{Biometrika}, 86(4), 948-955.
#'
#' @seealso \code{\link{fill_NA}} \code{\link{fill_NA_N}}
#'
#' @examples
#' library(miceFast)
#' set.seed(123)
#' data(air_miss)
#'
#' # Step 1: Generate m = 5 completed datasets using fill_NA with a stochastic model
#' completed <- lapply(1:5, function(i) {
#'   dat <- air_miss
#'   dat$Ozone <- fill_NA(
#'     x = dat,
#'     model = "lm_bayes",
#'     posit_y = "Ozone",
#'     posit_x = c("Solar.R", "Wind", "Temp")
#'   )
#'   dat
#' })
#'
#' # Step 2: Fit a model on each completed dataset
#' fits <- lapply(completed, function(d) {
#'   lm(Ozone ~ Solar.R + Wind + Temp, data = d)
#' })
#'
#' # Step 3: Pool using Rubin's rules
#' pool(fits)
#'
#' @export
pool <- function(fits, dfcom = NULL) {
  stopifnot(is.list(fits))
  if (length(fits) < 2L) {
    stop("At least 2 fitted models are required (length(fits) >= 2).", call. = FALSE)
  }

  m <- length(fits)

  # Extract coefficients and vcov from each fit
  coefs <- tryCatch(
    lapply(fits, stats::coef),
    error = function(e) {
      stop("All models must support coef(). Error: ", e$message, call. = FALSE)
    }
  )
  vcovs <- tryCatch(
    lapply(fits, stats::vcov),
    error = function(e) {
      stop("All models must support vcov(). Error: ", e$message, call. = FALSE)
    }
  )

  # Parameter names and dimension

  p <- length(coefs[[1]])
  nms <- names(coefs[[1]])
  if (is.null(nms)) nms <- paste0("V", seq_len(p))

  if (!all(vapply(coefs, length, integer(1)) == p)) {
    stop("All models must have the same number of coefficients.", call. = FALSE)
  }

  # Stack coefficients into m x p matrix
  Q <- do.call(rbind, coefs)

  # Pooled estimate: Qbar (Rubin 3.1.2)
  qbar <- colMeans(Q)

  # Within-imputation variance: Ubar (Rubin 3.1.3)
  ubar_mat <- Reduce("+", vcovs) / m

  # Between-imputation variance: B (Rubin 3.1.4)
  Q_centered <- sweep(Q, 2, qbar)
  b_mat <- crossprod(Q_centered) / (m - 1)

  # Total variance: T = Ubar + (1 + 1/m) * B (Rubin 3.1.5)
  total_mat <- ubar_mat + (1 + 1 / m) * b_mat

  # Diagonal elements for scalar per-parameter summaries
  ubar_diag <- diag(ubar_mat)
  b_diag <- diag(b_mat)
  t_diag <- diag(total_mat)

  # Relative increase in variance (Rubin 3.1.7)
  riv <- (1 + 1 / m) * b_diag / ubar_diag

  # Proportion of total variance attributable to missingness
  lambda <- (1 + 1 / m) * b_diag / t_diag

  # Complete-data degrees of freedom
  if (is.null(dfcom)) {
    dfcom_val <- tryCatch(
      {
        dfs <- vapply(fits, function(f) {
          r <- stats::df.residual(f)
          if (is.null(r)) NA_real_ else as.double(r)
        }, double(1))
        if (any(is.na(dfs))) Inf else mean(dfs)
      },
      error = function(e) Inf
    )
  } else {
    dfcom_val <- as.double(dfcom)
    stopifnot(dfcom_val > 0)
  }

  # Degrees of freedom: Barnard-Rubin (1999)
  df_old <- (m - 1) / lambda^2

  if (is.finite(dfcom_val)) {
    df_obs <- ((dfcom_val + 1) / (dfcom_val + 3)) * dfcom_val * (1 - lambda)
    df <- (df_old * df_obs) / (df_old + df_obs)
  } else {
    df <- df_old
  }

  # Standard errors, t-statistics, p-values
  se <- sqrt(t_diag)
  statistic <- qbar / se
  p_value <- 2 * stats::pt(-abs(statistic), df = df)

  # Fraction of missing information (Rubin 3.1.10)
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)

  # 95% confidence intervals
  crit <- stats::qt(0.975, df = df)
  conf_low <- qbar - crit * se
  conf_high <- qbar + crit * se

  result <- data.frame(
    term = nms,
    m = m,
    estimate = qbar,
    std.error = se,
    statistic = statistic,
    p.value = p_value,
    df = df,
    riv = riv,
    lambda = lambda,
    fmi = fmi,
    ubar = ubar_diag,
    b = b_diag,
    t = t_diag,
    dfcom = dfcom_val,
    conf.low = conf_low,
    conf.high = conf_high,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(result) <- c("miceFast_pool", "data.frame")
  result
}

#' Print method for pooled MI results
#'
#' @description Prints a concise summary of pooled multiple imputation results,
#'   showing the key inference columns: estimate, std.error, statistic, df, and p.value.
#' @param x an object of class \code{"miceFast_pool"}, as returned by \code{\link{pool}}.
#' @param ... additional arguments (currently ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.miceFast_pool <- function(x, ...) {
  cat("Pooled results from", x$m[1], "imputed datasets\n")
  cat("Rubin's rules with Barnard-Rubin df adjustment\n\n")
  display <- data.frame(
    term = x$term,
    estimate = x$estimate,
    std.error = x$std.error,
    statistic = x$statistic,
    df = x$df,
    p.value = x$p.value,
    stringsAsFactors = FALSE
  )
  print(display, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Summary method for pooled MI results
#'
#' @description Displays the full pooling diagnostics including within- and
#'   between-imputation variance, relative increase in variance, proportion of
#'   variance due to missingness, fraction of missing information, and
#'   confidence intervals.
#' @param object an object of class \code{"miceFast_pool"}, as returned by \code{\link{pool}}.
#' @param ... additional arguments (currently ignored).
#' @return Invisibly returns the full diagnostics data.frame.
#' @export
summary.miceFast_pool <- function(object, ...) {
  cat("Pooled results from", object$m[1], "imputed datasets\n")
  cat("Rubin's rules with Barnard-Rubin df adjustment\n")
  cat("Complete-data df:", object$dfcom[1], "\n\n")
  cat("Coefficients:\n")
  coef_df <- data.frame(
    term = object$term,
    estimate = object$estimate,
    std.error = object$std.error,
    statistic = object$statistic,
    df = object$df,
    p.value = object$p.value,
    `conf.low` = object$conf.low,
    `conf.high` = object$conf.high,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  print(coef_df, row.names = FALSE, digits = 4)
  cat("\nPooling diagnostics:\n")
  diag_df <- data.frame(
    term = object$term,
    ubar = object$ubar,
    b = object$b,
    t = object$t,
    riv = object$riv,
    lambda = object$lambda,
    fmi = object$fmi,
    stringsAsFactors = FALSE
  )
  print(diag_df, row.names = FALSE, digits = 4)
  invisible(object)
}
