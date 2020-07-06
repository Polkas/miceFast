
#' \code{fill_NA_N} function for the multiple imputations purpose.
#'
#' @description
#' Multiple imputations to fill the missing data.
#' Non missing independent variables are used to approximate a missing observations for a dependent variable.
#' Quantitative models were built under Rcpp packages and the C++ library Armadillo.
#'
#' @param x a numeric matrix or data.frame/data.table (factor/character/numeric/logical) - variables
#' @param model a character - posibble options ("lm_bayes","lm_noise","pmm")
#' @param posit_y an integer/character - a position/name of dependent variable
#' @param posit_x an integer/character vector - positions/names of independent variables
#' @param w  a numeric vector - a weighting variable - only positive values, Default: NULL
#' @param k an integer - a number of multiple imputations or for pmm a number of closest points from which a one random value is taken, Default:10
#' @param times deprecated
#' @param logreg a boolean - if dependent variable has log-normal distribution (numeric). If TRUE log-regression is evaluated and then returned exponential of results., Default: FALSE
#' @param ridge a numeric - a value added to diagonal elements of the x'x matrix, Default:1e-5
#'
#' @return load imputations in a numeric/character/factor (similar to the input type) vector format
#'
#' @note
#' There is assumed that users add the intercept by their own.
#' The miceFast module provides the most efficient environment, the second recommended option is to use data.table and the numeric matrix data type.
#' The lda model is assessed only if there are more than 15 complete observations
#' and for the lms models if number of variables is smaller than number of observations.
#'
#' @seealso \code{\link{fill_NA}} \code{\link{VIF}}
#'
#' @importFrom lifecycle deprecate_warn
#'
#' @examples
#' library(miceFast)
#' library(dplyr)
#' library(data.table)
#' ### Data
#' # airquality dataset with additional variables
#' data(air_miss)
#'
#' ### Intro: data.table
#' # IMPUTATIONS
#' # Imputations with a grouping option (models are separately assessed for each group)
#' # taking into account provided weights
#' air_miss[, Solar_R_imp := fill_NA_N(
#'   x = .SD,
#'   model = "lm_bayes",
#'   posit_y = "Solar.R",
#'   posit_x = c("Wind", "Temp", "Intercept"),
#'   w = .SD[["weights"]],
#'   times = 100
#' ), by = .(groups)] %>%
#'   # Imputations - discrete variable
#'   .[, x_character_imp := fill_NA(
#'     x = .SD,
#'     model = "lda",
#'     posit_y = "x_character",
#'     posit_x = c("Wind", "Temp", "groups")
#'   )] %>%
#'   # logreg was used because almost log-normal distribution of Ozone
#'   # imputations around mean
#'   .[, Ozone_imp1 := fill_NA(
#'     x = .SD,
#'     model = "lm_bayes",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept"),
#'     logreg = TRUE
#'   )] %>%
#'   # imputations using positions - Intercept, Temp
#'   .[, Ozone_imp2 := fill_NA(
#'     x = .SD,
#'     model = "lm_bayes",
#'     posit_y = 1,
#'     posit_x = c(4, 6),
#'     logreg = TRUE
#'   )] %>%
#'   # model with a factor independent variable
#'   # multiple imputations (average of x30 imputations)
#'   # with a factor independent variable, weights and logreg options
#'   .[, Ozone_imp3 := fill_NA_N(
#'     x = .SD,
#'     model = "lm_noise",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .SD[["weights"]],
#'     logreg = TRUE,
#'     times = 30
#'   )] %>%
#'   .[, Ozone_imp4 := fill_NA_N(
#'     x = .SD,
#'     model = "lm_bayes",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .SD[["weights"]],
#'     logreg = TRUE,
#'     times = 30
#'   )] %>%
#'   .[, Ozone_imp5 := fill_NA(
#'     x = .SD,
#'     model = "lm_pred",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .SD[["weights"]],
#'     logreg = TRUE
#'   ), .(groups)] %>%
#'   .[, Ozone_imp6 := fill_NA_N(
#'     x = .SD,
#'     model = "pmm",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .SD[["weights"]],
#'     logreg = TRUE,
#'     k = 10
#'   ), .(groups)] %>%
#'
#'   # Average of a few methods
#'   .[, Ozone_imp_mix := apply(.SD, 1, mean), .SDcols = Ozone_imp1:Ozone_imp6] %>%
#'
#'   # Protecting against collinearity or low number of observations - across small groups
#'   # Be carful when using a data.table grouping option
#'   # because of lack of protection against collinearity or low number of observations.
#'   # There could be used a tryCatch(fill_NA(...),error=function(e) return(...))
#'
#'   .[, Ozone_chac_imp := tryCatch(fill_NA(
#'     x = .SD,
#'     model = "lda",
#'     posit_y = "Ozone_chac",
#'     posit_x = c(
#'       "Intercept",
#'       "Month",
#'       "Day",
#'       "Temp",
#'       "x_character_imp"
#'     ),
#'     w = .SD[["weights"]]
#'   ),
#'   error = function(e) .SD[["Ozone_chac"]]
#'   ), .(groups)]
#'
#' # Sample of results
#' air_miss[which(is.na(air_miss[, 1]))[1:5], ]
#'
#' ### Intro: dplyr
#' # IMPUTATIONS
#' air_miss <- air_miss %>%
#'   # Imputations with a grouping option (models are separately assessed for each group)
#'   # taking into account provided weights
#'   group_by(groups) %>%
#'   do(mutate(., Solar_R_imp = fill_NA(
#'     x = .,
#'     model = "lm_pred",
#'     posit_y = "Solar.R",
#'     posit_x = c("Wind", "Temp", "Intercept"),
#'     w = .[["weights"]]
#'   ))) %>%
#'   ungroup() %>%
#'   # Imputations - discrete variable
#'   mutate(x_character_imp = fill_NA(
#'     x = .,
#'     model = "lda",
#'     posit_y = "x_character",
#'     posit_x = c("Wind", "Temp")
#'   )) %>%
#'   # logreg was used because almost log-normal distribution of Ozone
#'   # imputations around mean
#'   mutate(Ozone_imp1 = fill_NA(
#'     x = .,
#'     model = "lm_bayes",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept"),
#'     logreg = TRUE
#'   )) %>%
#'   # imputations using positions - Intercept, Temp
#'   mutate(Ozone_imp2 = fill_NA(
#'     x = .,
#'     model = "lm_bayes",
#'     posit_y = 1,
#'     posit_x = c(4, 6),
#'     logreg = TRUE
#'   )) %>%
#'   # multiple imputations (average of x30 imputations)
#'   # with a factor independent variable, weights and logreg options
#'   mutate(Ozone_imp3 = fill_NA_N(
#'     x = .,
#'     model = "lm_noise",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .[["weights"]],
#'     logreg = TRUE,
#'     times = 30
#'   )) %>%
#'   mutate(Ozone_imp4 = fill_NA_N(
#'     x = .,
#'     model = "lm_bayes",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .[["weights"]],
#'     logreg = TRUE,
#'     times = 30
#'   )) %>%
#'   group_by(groups) %>%
#'   do(mutate(., Ozone_imp5 = fill_NA(
#'     x = .,
#'     model = "lm_pred",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .[["weights"]],
#'     logreg = TRUE
#'   ))) %>%
#'   do(mutate(., Ozone_imp6 = fill_NA_N(
#'     x = .,
#'     model = "pmm",
#'     posit_y = "Ozone",
#'     posit_x = c("Intercept", "x_character_imp", "Wind", "Temp"),
#'     w = .[["weights"]],
#'     logreg = TRUE,
#'     k = 20
#'   ))) %>%
#'   ungroup() %>%
#'   # Average of a few methods
#'   mutate(Ozone_imp_mix = rowMeans(select(., starts_with("Ozone_imp")))) %>%
#'
#'   # Protecting against collinearity or low number of observations - across small groups
#'   # Be carful when using a grouping option
#'   # because of lack of protection against collinearity or low number of observations.
#'   # There could be used a tryCatch(fill_NA(...),error=function(e) return(...))
#'   group_by(groups) %>%
#'   do(mutate(., Ozone_chac_imp = tryCatch(fill_NA(
#'     x = .,
#'     model = "lda",
#'     posit_y = "Ozone_chac",
#'     posit_x = c(
#'       "Intercept",
#'       "Month",
#'       "Day",
#'       "Temp",
#'       "x_character_imp"
#'     ),
#'     w = .[["weights"]]
#'   ),
#'   error = function(e) .[["Ozone_chac"]]
#'   ))) %>%
#'   ungroup()
#'
#' # Sample of results
#' air_miss[which(is.na(air_miss[, 1]))[1:5], ]
#' @name fill_NA_N
#'
#' @export

fill_NA_N <- function(x, model, posit_y, posit_x, w = NULL, logreg = FALSE, k = 10, ridge = 1e-6, times = deprecated()) {
  if (inherits(x, "data.frame") || inherits(x, "matrix") || inherits(x, "data.table")) {
    if (is_present(times)) {
      deprecate_warn("0.6.0", "miceFast::fill_NA_N(times=)", "miceFast::fill_NA_N(k=)")
      k <- times
    }
    UseMethod("fill_NA_N")
  } else {
    stop("wrong data type - it should be data.frame, matrix or data.table")
  }
}

#' @describeIn fill_NA_N s3 method for data.frame

fill_NA_N.data.frame <- function(x, model, posit_y, posit_x, w = NULL, logreg = FALSE, k = 10, ridge = 1e-6, times = deprecated()) {
  ww <- if (is.null(w)) vector() else w

  if (posit_y %in% posit_x) {
    stop("the same variable is dependent and indepentent")
  }
  model <- match.arg(model, c("lm_bayes", "lm_noise", "pmm"))


  cols <- colnames(x)

  if (is.character(posit_x)) {
    posit_x <- pmatch(posit_x, cols)
    posit_x <- posit_x[!is.na(posit_x)]
    if (length(posit_x) == 0) stop("posit_x is empty")
  }

  if (is.character(posit_y)) {
    posit_y <- pmatch(posit_y, cols)
    if (length(posit_y) == 0) stop("posit_y is empty")
  }

  yy <- x[[posit_y]]

  yy_class <- class(yy)

  is_factor_y <- yy_class == "factor"
  is_character_y <- yy_class == "character"
  is_numeric_y <- (yy_class == "numeric") || (yy_class == "integer") || (yy_class == "logical") || (yy_class == "logical")

  all_pos_y <- FALSE
  if (is_numeric_y) {
    all_pos_y <- !any(yy < 0, na.rm = TRUE)
  }

  if ((is_character_y || is_factor_y || (model == "lda")) && logreg) {
    stop("logreg works only for a non-negative numeric dependent variable and lm models")
  } else if (all_pos_y && logreg) {
    yy <- log(yy + 1e-8)
  }

  x_small <- x[, posit_x]
  types <- lapply(x_small, class)
  x_ncols <- length(posit_x)
  p_x_factor_character <- which(unlist(lapply(types, function(i) !all(is.na(match(c("factor", "character"), i))))))
  len_p_x_factor_character <- length(p_x_factor_character)

  xx <- vector("list", 2)

  if (len_p_x_factor_character > 0) {
    posit_fc <- posit_x[p_x_factor_character]
    x_fc <- x[, posit_fc]
    x_fc <- model.matrix.lm(~., x_fc, na.action = "na.pass")[, -1]
    xx[[1]] <- x_fc
  }

  if (x_ncols > len_p_x_factor_character) {
    posit_ni <- setdiff(posit_x, posit_x[p_x_factor_character])
    x_ni <- as.matrix(x[, posit_ni])
    xx[[2]] <- x_ni
  }

  xx <- do.call(cbind, xx[!is.null(xx)])

  if (is_factor_y) {
    l <- levels(yy)
    yy <- as.numeric(yy)
    f <- round(fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge))
    f[f <= 0] <- 1
    f[f > length(l)] <- length(l)
    ff <- factor(l[f])
  } else if (is_character_y) {
    yy <- factor(yy)
    l <- levels(yy)
    yy <- as.numeric(yy)
    yy <- yy
    f <- round(fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge))
    f[f <= 0] <- 1
    f[f > length(l)] <- length(l)
    ff <- l[f]
  } else if (is_numeric_y) {
    yy <- as.numeric(yy)
    ff <- fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge)
    if (logreg && (model != "lda")) {
      ff <- exp(ff)
    }
  }

  attr(ff, "dim") <- attributes(ff)$dim[1]

  return(ff)
}

#' @describeIn fill_NA_N S3 method for data.table

fill_NA_N.data.table <- function(x, model, posit_y, posit_x, w = NULL, logreg = FALSE, k = 10, ridge = 1e-6, times = deprecated()) {
  ww <- if (is.null(w)) vector() else w
  if (posit_y %in% posit_x) {
    stop("the same variable is dependent and indepentent")
  }
  model <- match.arg(model, c("lm_bayes", "lm_noise", "pmm"))

  cols <- colnames(x)

  if (is.character(posit_x)) {
    posit_x <- pmatch(posit_x, cols)
    posit_x <- posit_x[!is.na(posit_x)]
    if (length(posit_x) == 0) stop("posit_x is empty")
  }

  if (is.character(posit_y)) {
    posit_y <- pmatch(posit_y, cols)
    if (length(posit_y) == 0) stop("posit_y is empty")
  }

  yy <- x[[posit_y]]

  yy_class <- class(yy)

  is_factor_y <- yy_class == "factor"
  is_character_y <- yy_class == "character"
  is_numeric_y <- (yy_class == "numeric") || (yy_class == "integer") || (yy_class == "logical") || (yy_class == "logical")

  all_pos_y <- FALSE
  if (is_numeric_y) {
    all_pos_y <- !any(yy < 0, na.rm = TRUE)
  }

  if ((is_character_y || is_factor_y || (model == "lda")) && logreg) {
    stop("logreg works only for a non-negative numeric dependent variable and lm models")
  } else if (all_pos_y && logreg) {
    yy <- log(yy + 1e-8)
  }

  x_small <- x[, posit_x, with = FALSE]
  types <- lapply(x_small, class)
  x_ncols <- length(posit_x)
  p_x_factor_character <- which(unlist(lapply(types, function(i) !all(is.na(match(c("factor", "character"), i))))))
  len_p_x_factor_character <- length(p_x_factor_character)

  xx <- vector("list", 2)

  if (len_p_x_factor_character > 0) {
    posit_fc <- posit_x[p_x_factor_character]
    x_fc <- x[, posit_fc, with = FALSE]
    x_fc <- model.matrix.lm(~., x_fc, na.action = "na.pass")[, -1]
    xx[[1]] <- x_fc
  }

  if (x_ncols > len_p_x_factor_character) {
    posit_ni <- setdiff(posit_x, posit_x[p_x_factor_character])
    x_ni <- as.matrix(x[, posit_ni, with = FALSE])
    xx[[2]] <- x_ni
  }

  xx <- do.call(cbind, xx[!is.null(xx)])

  if (is_factor_y) {
    l <- levels(yy)
    yy <- as.numeric(yy)
    f <- round(fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge))
    f[f <= 0] <- 1
    f[f > length(l)] <- length(l)
    ff <- factor(l[f])
  } else if (is_character_y) {
    yy <- factor(yy)
    l <- levels(yy)
    yy <- as.numeric(yy)
    yy <- yy
    f <- round(fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge))
    f[f <= 0] <- 1
    f[f > length(l)] <- length(l)
    ff <- l[f]
  } else if (is_numeric_y) {
    yy <- as.numeric(yy)
    ff <- fill_NA_N_(cbind(yy, xx), model, 1, 2:(ncol(xx) + 1), ww, k, ridge)
    if (logreg && (model != "lda")) {
      ff <- exp(ff)
    }
  }

  attr(ff, "dim") <- attributes(ff)$dim[1]

  return(ff)
}

#' @describeIn fill_NA_N S3 method for matrix

fill_NA_N.matrix <- function(x, model, posit_y, posit_x, w = NULL, logreg = FALSE, k = 10, ridge = 1e-6, times = deprecated()) {
  ww <- if (is.null(w)) vector() else w
  if (posit_y %in% posit_x) {
    stop("the same variable is dependent and indepentent")
  }
  model <- match.arg(model, c("lm_bayes", "lm_noise", "pmm"))

  cols <- colnames(x)

  if (is.character(posit_x)) {
    posit_x <- pmatch(posit_x, cols)
    posit_x <- posit_x[!is.na(posit_x)]
    if (length(posit_x) == 0) stop("posit_x is empty")
  }

  if (is.character(posit_y)) {
    posit_y <- pmatch(posit_y, cols)
    if (length(posit_y) == 0) stop("posit_y is empty")
  }

  all_pos_y <- !any(x[[posit_y]] < 0, na.rm = TRUE)
  logreg_con <- logreg && all_pos_y && (model != "lda")

  if (logreg_con) {
    x[[posit_y]] <- log(x[[posit_y]] + 1e-8)
  }
  ff <- fill_NA_N_(x, model, posit_y, posit_x, ww, k, ridge)
  if (logreg_con) {
    ff <- exp(ff)
  }

  attr(ff, "dim") <- attributes(ff)$dim[1]

  return(ff)
}
