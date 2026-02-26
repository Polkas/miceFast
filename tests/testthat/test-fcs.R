# Tests for iterative FCS (chained equations) pattern
# Compares miceFast FCS loop with mice on the air_miss dataset

library(data.table)

# --- Helper: data.table FCS function -----------------------------------------

fcs_dt <- function(dt, n_iter = 5) {
  d <- copy(dt)
  na_ozone <- is.na(d$Ozone)
  na_solar <- is.na(d[["Solar.R"]])
  d <- naive_fill_NA(d)
  for (iter in seq_len(n_iter)) {
    set(d, which(na_ozone), "Ozone", NA_real_)
    d[, Ozone := fill_NA(.SD, "lm_bayes", "Ozone", c("Solar.R", "Wind", "Temp"))]
    set(d, which(na_solar), "Solar.R", NA_real_)
    d[, Solar.R := fill_NA(.SD, "lm_bayes", "Solar.R", c("Ozone", "Wind", "Temp", "Intercept"))]
  }
  d
}

# --- Helper: data.frame FCS function -----------------------------------------

fcs_df <- function(df, n_iter = 5) {
  na_ozone <- which(is.na(df$Ozone))
  na_solar <- which(is.na(df$Solar.R))
  d <- naive_fill_NA(df)
  for (iter in seq_len(n_iter)) {
    d$Ozone[na_ozone] <- NA
    d$Ozone <- fill_NA(d, "lm_bayes", "Ozone", c("Solar.R", "Wind", "Temp"))
    d$Solar.R[na_solar] <- NA
    d$Solar.R <- fill_NA(d, "lm_bayes", "Solar.R", c("Ozone", "Wind", "Temp", "Intercept"))
  }
  d
}

# --- Helper: OOP FCS function ------------------------------------------------

fcs_oop <- function(mat, n_iter = 5) {
  na_ozone <- is.na(mat[, 1])
  na_solar <- is.na(mat[, 2])
  mat_init <- mat + 0
  for (j in seq_len(ncol(mat_init))) {
    na_j <- is.na(mat_init[, j])
    if (any(na_j)) mat_init[na_j, j] <- mean(mat_init[!na_j, j])
  }
  m <- new(miceFast)
  m$set_data(mat_init)
  for (iter in seq_len(n_iter)) {
    col1 <- m$get_data()[, 1]; col1[na_ozone] <- NaN
    m$update_var(1, col1)
    m$update_var(1, m$impute("lm_bayes", 1, c(2, 3, 4))$imputations)
    col2 <- m$get_data()[, 2]; col2[na_solar] <- NaN
    m$update_var(2, col2)
    m$update_var(2, m$impute("lm_bayes", 2, c(1, 3, 4, 5))$imputations)
  }
  as.data.frame(m$get_data())
}

# =============================================================================
# Tests: Basic FCS correctness
# =============================================================================

test_that("FCS data.table: all NAs are imputed", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  set.seed(101)
  result <- fcs_dt(dt)
  expect_equal(sum(is.na(result$Ozone)), 0)
  expect_equal(sum(is.na(result$Solar.R)), 0)
})

test_that("FCS data.frame: all NAs are imputed", {
  data(air_miss, package = "miceFast")
  set.seed(102)
  result <- fcs_df(air_miss)
  expect_equal(sum(is.na(result$Ozone)), 0)
  expect_equal(sum(is.na(result$Solar.R)), 0)
})

test_that("FCS OOP matrix: all NaN are imputed", {
  data(air_miss, package = "miceFast")
  mat <- cbind(
    as.matrix(air_miss[, c("Ozone", "Solar.R", "Wind", "Temp")]),
    Intercept = 1
  )
  set.seed(103)
  result <- fcs_oop(mat)
  expect_equal(sum(is.nan(result[[1]])), 0)
  expect_equal(sum(is.nan(result[[2]])), 0)
})

test_that("FCS handles joint missingness (rows where both Ozone and Solar.R are NA)", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  joint_na <- which(is.na(dt$Ozone) & is.na(dt$Solar.R))
  expect_true(length(joint_na) > 0, info = "air_miss should have jointly missing rows")

  set.seed(104)
  result <- fcs_dt(dt)
  expect_false(anyNA(result$Ozone[joint_na]))
  expect_false(anyNA(result$Solar.R[joint_na]))
})

test_that("FCS does not modify the original data", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  dt_orig <- copy(dt)
  set.seed(105)
  fcs_dt(dt)
  expect_identical(dt, dt_orig)
})

test_that("FCS preserves non-missing values", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  obs_ozone <- which(!is.na(dt$Ozone))
  obs_solar <- which(!is.na(dt$Solar.R))
  orig_ozone <- dt$Ozone[obs_ozone]
  orig_solar <- dt$Solar.R[obs_solar]

  set.seed(106)
  result <- fcs_dt(dt)
  expect_equal(result$Ozone[obs_ozone], orig_ozone)
  expect_equal(result$Solar.R[obs_solar], orig_solar)
})

test_that("FCS n_iter = 1 still produces complete data", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  set.seed(107)
  result <- fcs_dt(dt, n_iter = 1)
  expect_equal(sum(is.na(result$Ozone)), 0)
  expect_equal(sum(is.na(result$Solar.R)), 0)
})

test_that("FCS produces different results across stochastic runs (lm_bayes)", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)
  set.seed(201)
  r1 <- fcs_dt(dt)
  set.seed(202)
  r2 <- fcs_dt(dt)
  # Imputed values should differ between seeds
  na_idx <- which(is.na(dt$Ozone))
  expect_false(identical(r1$Ozone[na_idx], r2$Ozone[na_idx]))
})

# =============================================================================
# Tests: MI + pool workflow
# =============================================================================

test_that("FCS MI workflow produces valid pooled results", {
  data(air_miss, package = "miceFast")
  dt <- as.data.table(air_miss)

  set.seed(301)
  completed <- lapply(1:10, function(i) fcs_dt(dt))
  fits <- lapply(completed, function(d) lm(Ozone ~ Wind + Temp, data = d))
  res <- pool(fits)

  expect_s3_class(res, "miceFast_pool")
  s <- summary(res)

  # Coefficients should be reasonable for Ozone ~ Wind + Temp
  # Wind should have a negative effect, Temp a positive effect
  expect_true(s[s$term == "Wind", "estimate"] < 0)
  expect_true(s[s$term == "Temp", "estimate"] > 0)

  # Standard errors should be positive and finite
  expect_true(all(s$std.error > 0))
  expect_true(all(is.finite(s$std.error)))

  # p-values should exist and be in [0, 1]
  expect_true(all(s$p.value >= 0 & s$p.value <= 1))
})

# =============================================================================
# Tests: Comparison with mice
# =============================================================================

test_that("FCS pooled estimates are similar to mice", {
  skip_if_not_installed("mice")

  data(air_miss, package = "miceFast")

  # --- mice ---
  mice_data <- air_miss[, c("Ozone", "Solar.R", "Wind", "Temp")]
  set.seed(42)
  mice_imp <- mice::mice(
    mice_data,
    m = 20,
    method = "norm",
    maxit = 5,
    printFlag = FALSE
  )
  mice_fits <- with(mice_imp, lm(Ozone ~ Wind + Temp))
  mice_pool <- mice::pool(mice_fits)
  mice_coefs <- summary(mice_pool)

  # --- miceFast FCS ---
  dt <- as.data.table(air_miss)
  set.seed(42)
  completed <- lapply(1:20, function(i) fcs_dt(dt))
  fits <- lapply(completed, function(d) lm(Ozone ~ Wind + Temp, data = d))
  mf_pool <- pool(fits)
  mf_coefs <- summary(mf_pool)

  # Point estimates should be in the same ballpark (within 50% or 15 units)
  for (trm in c("(Intercept)", "Wind", "Temp")) {
    mice_est <- mice_coefs[mice_coefs$term == trm, "estimate"]
    mf_est <- mf_coefs[mf_coefs$term == trm, "estimate"]
    expect_true(
      abs(mf_est - mice_est) < max(abs(mice_est) * 0.5, 15),
      info = sprintf(
        "%s: miceFast=%.3f vs mice=%.3f (diff=%.3f)",
        trm, mf_est, mice_est, mf_est - mice_est
      )
    )
  }

  # Sign of effects should match
  expect_equal(
    sign(mf_coefs[mf_coefs$term == "Wind", "estimate"]),
    sign(mice_coefs[mice_coefs$term == "Wind", "estimate"])
  )
  expect_equal(
    sign(mf_coefs[mf_coefs$term == "Temp", "estimate"]),
    sign(mice_coefs[mice_coefs$term == "Temp", "estimate"])
  )

  # Standard errors should be in the same order of magnitude
  for (trm in c("(Intercept)", "Wind", "Temp")) {
    mice_se <- mice_coefs[mice_coefs$term == trm, "std.error"]
    mf_se <- mf_coefs[mf_coefs$term == trm, "std.error"]
    ratio <- mf_se / mice_se
    expect_true(
      ratio > 0.25 && ratio < 4,
      info = sprintf(
        "%s SE ratio: miceFast/mice = %.2f (mf=%.3f, mice=%.3f)",
        trm, ratio, mf_se, mice_se
      )
    )
  }
})

test_that("FCS imputed means are close to mice imputed means", {
  skip_if_not_installed("mice")

  data(air_miss, package = "miceFast")

  # --- mice ---
  mice_data <- air_miss[, c("Ozone", "Solar.R", "Wind", "Temp")]
  set.seed(99)
  mice_imp <- mice::mice(
    mice_data,
    m = 20,
    method = "norm",
    maxit = 5,
    printFlag = FALSE
  )
  mice_completed <- mice::complete(mice_imp, action = "all")
  mice_ozone_mean <- mean(sapply(mice_completed, function(d) mean(d$Ozone)))
  mice_solar_mean <- mean(sapply(mice_completed, function(d) mean(d$Solar.R)))

  # --- miceFast FCS ---
  dt <- as.data.table(air_miss)
  set.seed(99)
  completed <- lapply(1:20, function(i) fcs_dt(dt))
  mf_ozone_mean <- mean(sapply(completed, function(d) mean(d$Ozone)))
  mf_solar_mean <- mean(sapply(completed, function(d) mean(d$Solar.R)))

  # Means should be within 20% of each other
  expect_true(
    abs(mf_ozone_mean - mice_ozone_mean) / mice_ozone_mean < 0.2,
    info = sprintf("Ozone mean: miceFast=%.2f vs mice=%.2f", mf_ozone_mean, mice_ozone_mean)
  )
  expect_true(
    abs(mf_solar_mean - mice_solar_mean) / mice_solar_mean < 0.2,
    info = sprintf("Solar.R mean: miceFast=%.2f vs mice=%.2f", mf_solar_mean, mice_solar_mean)
  )
})
