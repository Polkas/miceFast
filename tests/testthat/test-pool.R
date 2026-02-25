# ---- Setup ----
set.seed(12345)

# Create a dataset with missing values for testing
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + 2 * x1 + 3 * x2 + rnorm(n)
dat <- data.frame(y = y, x1 = x1, x2 = x2)

# Introduce MCAR missingness
miss_idx <- sample(seq_len(n), size = 40)
dat$y[miss_idx] <- NA

m <- 10

# ---- Tests for pool() ----

test_that("pool returns correct structure", {
  completed <- lapply(1:m, function(i) {
    d <- dat
    d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)), mean(d$y, na.rm = TRUE), sd(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) lm(y ~ x1 + x2, data = d))

  result <- pool(fits)

  expect_s3_class(result, "data.frame")
  expect_s3_class(result, "miceFast_pool")
  expect_equal(nrow(result), 3L)
  expect_named(result, c(
    "term", "m", "estimate", "std.error", "statistic", "p.value",
    "df", "riv", "lambda", "fmi", "ubar", "b", "t", "dfcom",
    "conf.low", "conf.high"
  ))
  expect_equal(result$term, c("(Intercept)", "x1", "x2"))
  expect_true(all(result$m == m))
  expect_true(all(result$std.error > 0))
  expect_true(all(result$df > 0))
  expect_true(all(result$t > 0))
  expect_true(all(result$ubar > 0))
  expect_true(all(result$b >= 0))
  expect_true(all(result$riv >= 0))
  expect_true(all(result$lambda >= 0 & result$lambda <= 1))
  expect_true(all(result$fmi >= 0 & result$fmi <= 1))
  expect_true(all(result$conf.low < result$estimate))
  expect_true(all(result$conf.high > result$estimate))
})

test_that("pool errors on less than 2 fits", {
  d <- dat
  d$y[is.na(d$y)] <- 0
  fit1 <- lm(y ~ x1 + x2, data = d)

  expect_error(pool(list(fit1)), "At least 2 fitted models")
  expect_error(pool(list()), "At least 2 fitted models")
})

test_that("pool errors on mismatched coefficients", {
  d <- dat
  d$y[is.na(d$y)] <- 0
  fit1 <- lm(y ~ x1 + x2, data = d)
  fit2 <- lm(y ~ x1, data = d)

  expect_error(pool(list(fit1, fit2)), "same number of coefficients")
})

test_that("pool works with glm", {
  dat_bin <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  dat_bin$y[sample(n, 20)] <- NA

  completed <- lapply(1:5, function(i) {
    d <- dat_bin
    d$y[is.na(d$y)] <- rbinom(sum(is.na(d$y)), 1, mean(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) glm(y ~ x1 + x2, data = d, family = binomial))

  result <- pool(fits)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_true(all(is.finite(result$estimate)))
  expect_true(all(is.finite(result$std.error)))
})

test_that("pool respects dfcom argument", {
  completed <- lapply(1:m, function(i) {
    d <- dat
    d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)), mean(d$y, na.rm = TRUE), sd(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) lm(y ~ x1 + x2, data = d))

  result_auto <- pool(fits)
  result_inf <- pool(fits, dfcom = Inf)
  result_small <- pool(fits, dfcom = 10)

  # With smaller dfcom the adjusted df should be smaller
  expect_true(all(result_small$df <= result_auto$df + 1e-10))
  # dfcom = Inf should give df_old = (m-1)/lambda^2
  expect_true(all(result_inf$dfcom == Inf))
  expect_true(all(is.finite(result_inf$df)))
})

test_that("pool Rubin's rules match manual calculation", {
  # Generate deterministic completed datasets for reproducibility
  set.seed(999)
  completed <- lapply(1:5, function(i) {
    d <- dat
    d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)), mean(d$y, na.rm = TRUE), sd(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) lm(y ~ x1 + x2, data = d))

  result <- pool(fits)
  m_val <- 5

  # Manual calculation
  coefs_mat <- do.call(rbind, lapply(fits, coef))
  vcovs_list <- lapply(fits, vcov)

  qbar_manual <- colMeans(coefs_mat)
  ubar_manual <- Reduce("+", vcovs_list) / m_val
  Q_centered <- sweep(coefs_mat, 2, qbar_manual)
  b_manual <- crossprod(Q_centered) / (m_val - 1)
  t_manual <- diag(ubar_manual) + (1 + 1 / m_val) * diag(b_manual)

  expect_equal(result$estimate, unname(qbar_manual), tolerance = 1e-12)
  expect_equal(result$ubar, unname(diag(ubar_manual)), tolerance = 1e-12)
  expect_equal(result$b, unname(diag(b_manual)), tolerance = 1e-12)
  expect_equal(result$t, unname(t_manual), tolerance = 1e-12)
  expect_equal(result$std.error, sqrt(unname(t_manual)), tolerance = 1e-12)
})

test_that("pool matches mice::pool results", {
  skip_if_not_installed("mice")

  # Save reference to our pool before mice potentially masks it
  miceFast_pool <- pool

  set.seed(42)

  # Use mice::nhanes as a clean test case
  nhanes <- mice::nhanes

  # Run mice to get m imputed datasets
  imp <- mice::mice(nhanes, m = 5, method = "norm.predict", seed = 42, printFlag = FALSE)

  # Fit model on each imputed dataset via mice's with()
  mice_fit <- with(imp, lm(chl ~ age + bmi))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  # Extract the same completed datasets and fit with base lm
  completed_list <- lapply(1:5, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) lm(chl ~ age + bmi, data = d))

  # Pool with miceFast
  mf_pooled <- miceFast_pool(base_fits)

  # Compare pooled estimates
  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)

  # Compare within-imputation variance
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)

  # Compare between-imputation variance
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)

  # Compare total variance
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)

  # Compare standard errors
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)

  # Compare degrees of freedom
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)

  # Compare riv (relative increase in variance)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)

  # Compare lambda
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)

  # Compare fmi
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)

  # Compare t-statistics
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)

  # Compare p-values
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool for glm (logistic regression)", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(101)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 + 0.8 * x1 - 0.6 * x2))
  dat_logistic <- data.frame(y = factor(y), x1 = x1, x2 = x2)
  dat_logistic$y[sample(n, 30)] <- NA

  imp <- mice::mice(dat_logistic, m = 5, method = c("logreg", "norm", "norm"),
                    seed = 101, printFlag = FALSE)

  mice_fit <- with(imp, glm(y ~ x1 + x2, family = binomial))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  completed_list <- lapply(1:5, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) {
    glm(y ~ x1 + x2, data = d, family = binomial)
  })

  mf_pooled <- miceFast_pool(base_fits)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool for glm (Poisson regression)", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(202)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n, 2, 0.5)
  y <- rpois(n, exp(0.5 + 0.3 * x1 + 0.2 * x2))
  dat_pois <- data.frame(y = y, x1 = x1, x2 = x2)
  dat_pois$x1[sample(n, 25)] <- NA

  imp <- mice::mice(dat_pois, m = 7, method = c("pmm", "norm", "norm"),
                    seed = 202, printFlag = FALSE)

  mice_fit <- with(imp, glm(y ~ x1 + x2, family = poisson))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  completed_list <- lapply(1:7, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) {
    glm(y ~ x1 + x2, data = d, family = poisson)
  })

  mf_pooled <- miceFast_pool(base_fits)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool for lm with interactions", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(303)
  nhanes2 <- mice::nhanes2

  imp <- mice::mice(nhanes2, m = 5, method = "pmm", seed = 303, printFlag = FALSE)

  mice_fit <- with(imp, lm(chl ~ bmi * hyp))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  completed_list <- lapply(1:5, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) {
    lm(chl ~ bmi * hyp, data = d)
  })

  mf_pooled <- miceFast_pool(base_fits)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool with many imputations (m=20)", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(404)
  nhanes <- mice::nhanes

  imp <- mice::mice(nhanes, m = 20, method = "norm", seed = 404, printFlag = FALSE)

  mice_fit <- with(imp, lm(chl ~ age + bmi))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  completed_list <- lapply(1:20, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) {
    lm(chl ~ age + bmi, data = d)
  })

  mf_pooled <- miceFast_pool(base_fits)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool for weighted lm (with dfcom)", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(505)
  n <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 1.5 * x1 - 0.5 * x2 + rnorm(n)
  w <- runif(n, 0.5, 2)
  dat_wt <- data.frame(y = y, x1 = x1, x2 = x2, w = w)
  dat_wt$y[sample(n, 20)] <- NA

  imp <- mice::mice(dat_wt, m = 5, seed = 505, printFlag = FALSE)

  # mice with() does not propagate weights easily, so we extract datasets
  # and fit weighted lm manually on both sides
  completed_list <- lapply(1:5, function(i) mice::complete(imp, i))

  # Fit weighted lm from mice's completed datasets
  mice_fits_raw <- lapply(completed_list, function(d) {
    lm(y ~ x1 + x2, data = d, weights = w)
  })

  # Pool via mice (using as.mira to wrap base fits)
  mice_pooled <- mice::pool(mice::as.mira(mice_fits_raw))
  mice_summary <- summary(mice_pooled)

  # Pool via miceFast
  mf_pooled <- miceFast_pool(mice_fits_raw)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("pool matches mice::pool for multivariate missingness", {
  skip_if_not_installed("mice")

  miceFast_pool <- pool

  set.seed(606)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n, 5, 2)
  x3 <- rbinom(n, 1, 0.4)
  y <- 3 + x1 - 2 * x2 + 1.5 * x3 + rnorm(n)
  dat_multi <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
  # Introduce missingness in multiple variables
  dat_multi$y[sample(n, 40)] <- NA
  dat_multi$x1[sample(n, 30)] <- NA
  dat_multi$x2[sample(n, 25)] <- NA

  imp <- mice::mice(dat_multi, m = 10, seed = 606, printFlag = FALSE)

  mice_fit <- with(imp, lm(y ~ x1 + x2 + x3))
  mice_pooled <- mice::pool(mice_fit)
  mice_summary <- summary(mice_pooled)

  completed_list <- lapply(1:10, function(i) mice::complete(imp, i))
  base_fits <- lapply(completed_list, function(d) {
    lm(y ~ x1 + x2 + x3, data = d)
  })

  mf_pooled <- miceFast_pool(base_fits)

  expect_equal(mf_pooled$estimate, as.numeric(mice_summary$estimate), tolerance = 1e-8)
  expect_equal(mf_pooled$std.error, as.numeric(mice_summary$std.error), tolerance = 1e-8)
  expect_equal(mf_pooled$ubar, mice_pooled$pooled$ubar, tolerance = 1e-8)
  expect_equal(mf_pooled$b, mice_pooled$pooled$b, tolerance = 1e-8)
  expect_equal(mf_pooled$t, mice_pooled$pooled$t, tolerance = 1e-8)
  expect_equal(mf_pooled$df, mice_pooled$pooled$df, tolerance = 1e-6)
  expect_equal(mf_pooled$riv, mice_pooled$pooled$riv, tolerance = 1e-8)
  expect_equal(mf_pooled$lambda, mice_pooled$pooled$lambda, tolerance = 1e-8)
  expect_equal(mf_pooled$fmi, mice_pooled$pooled$fmi, tolerance = 1e-6)
  expect_equal(mf_pooled$statistic, as.numeric(mice_summary$statistic), tolerance = 1e-8)
  expect_equal(mf_pooled$p.value, as.numeric(mice_summary$p.value), tolerance = 1e-6)
})

test_that("print.mipo produces output", {
  completed <- lapply(1:5, function(i) {
    d <- dat
    d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)), mean(d$y, na.rm = TRUE), sd(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) lm(y ~ x1 + x2, data = d))
  result <- pool(fits)

  out <- capture.output(print(result))
  expect_true(any(grepl("Pooled results from 5 imputed datasets", out)))
  expect_true(any(grepl("Rubin's rules", out)))
  expect_true(any(grepl("\\(Intercept\\)", out)))

  # print returns invisibly
  expect_identical(print(result), result)
})

test_that("summary.mipo shows full diagnostics", {
  completed <- lapply(1:5, function(i) {
    d <- dat
    d$y[is.na(d$y)] <- rnorm(sum(is.na(d$y)), mean(d$y, na.rm = TRUE), sd(d$y, na.rm = TRUE))
    d
  })
  fits <- lapply(completed, function(d) lm(y ~ x1 + x2, data = d))
  result <- pool(fits)

  out <- capture.output(summary(result))
  expect_true(any(grepl("Pooled results from 5 imputed datasets", out)))
  expect_true(any(grepl("Complete-data df:", out)))
  expect_true(any(grepl("Coefficients:", out)))
  expect_true(any(grepl("Pooling diagnostics:", out)))
  expect_true(any(grepl("conf.low", out)))
  expect_true(any(grepl("fmi", out)))

  # summary returns invisibly
  expect_identical(summary(result), result)
})
