library(dplyr)
library(data.table)

set.seed(1234)

power <- 4 # power of 10 - number of observations - should be adjusted to a computer capabilities
nr_var <- 7 # CHANGE - only if you generate a bigger corr matrix:  number of variables - independent and one dependent
grs <- max(c(10**(power - 2), 10)) # grouping variable - number of groups

## generete example - data
## positive-defined correlation matrix
cors <- matrix(c(
  1, 0.76, 0.60, 0.44, 0.45, 0.45, 0.37,
  NA, 1, 0.25, 0.09, 0.13, 0.15, 0.17,
  NA, NA, 1, 0.17, 0.17, 0.12, 0.09,
  NA, NA, NA, 1, 0.13, 0.16, 0.11,
  NA, NA, NA, NA, 1, 0.16, 0.21,
  NA, NA, NA, NA, NA, 1, 0.16,
  NA, NA, NA, NA, NA, NA, 1
), 7, 7, byrow = T)
cors[lower.tri(cors)] <- t(cors)[lower.tri(cors)]
##

model <- new(corrData, 10, 10^power, rep(0, nr_var), cors)
data_bin <- model$fill("binom")
data_disc <- model$fill("discrete")
data_con <- model$fill("contin")

rm(model)

n_vars <- ncol(cors)

posit_y <- 1
posit_x <- 2:(n_vars - 2)
posit_w <- n_vars - 1
posit_grs <- n_vars
posit_NA <- n_vars + 1

## NA index

index_NA <- 1:nrow(data_con) %in% sample(1:nrow(data_con), 10^(power - 1))

fill_by_NA <- function(v, index_NA) {
  v[index_NA] <- NA
  v
}

######################

group_d <- floor(pnorm(data_disc[, posit_grs]) * grs)
group_c <- floor(pnorm(data_con[, posit_grs]) * grs)
group_b <- floor(pnorm(data_bin[, posit_grs]) * grs)

w_d <- abs(data_disc[, posit_w])
w_c <- abs(data_con[, posit_w])
w_b <- abs(data_bin[, posit_w])

data_disc_NA <- cbind(fill_by_NA(data_disc[, posit_y], index_NA), data_disc[, posit_x], w_d, group_d, index_NA)
data_con_NA <- cbind(fill_by_NA(data_con[, posit_y], index_NA), data_con[, posit_x], w_c, group_c, index_NA)
data_bin_NA <- cbind(fill_by_NA(data_bin[, posit_y], index_NA), data_bin[, posit_x], w_b, group_b, index_NA)

colnames(data_bin_NA) <- c("y", paste0("x", posit_x), "weights", "group", "index_NA")
colnames(data_disc_NA) <- c("y", paste0("x", posit_x), "weights", "group", "index_NA")
colnames(data_con_NA) <- c("y", paste0("x", posit_x), "weights", "group", "index_NA")


#################################

testthat::test_that("fill_NA lm_bayes accuracy - matrix", {
  pred_lm_bayes <- fill_NA(
    x = as.matrix(data_con_NA),
    model = "lm_bayes",
    posit_y = 1,
    posit_x = c(2, 3, 4)
  )
  testthat::expect_true(cor(cbind(pred_lm_bayes[index_NA], data_con[index_NA, 1]))[1, 2] > 0.5)
})

testthat::test_that("fill_NA_N lm_bayes accuracy", {
  pred_lm_bayes_N <- fill_NA_N(
    x = as.matrix(data_con_NA),
    model = "lm_bayes",
    posit_y = 1,
    posit_x = c(2, 3, 4)
  )
  testthat::expect_true(cor(cbind(pred_lm_bayes_N[index_NA], data_con[index_NA, 1]))[1, 2] > 0.5)
})

testthat::test_that("fill_NA lm_noise accuracy", {
  pred_lm_noise <- fill_NA(
    x = as.matrix(data_con_NA),
    model = "lm_noise",
    posit_y = 1,
    posit_x = c(2, 3, 4)
  )
  testthat::expect_true(cor(cbind(pred_lm_noise[index_NA], data_con[index_NA, 1]))[1, 2] > 0.5)
})

testthat::test_that("fill_NA_N lm_noise accuracy", {
  pred_lm_noise_N <- fill_NA_N(
    x = as.matrix(data_con_NA),
    model = "lm_noise",
    posit_y = 1,
    posit_x = c(2, 3, 4),
    k = 10
  )
  testthat::expect_true(cor(cbind(pred_lm_noise_N[index_NA], data_con[index_NA, 1]))[1, 2] > 0.5)
})

testthat::test_that("fill_NA lda accuracy", {
  pred_lda <- fill_NA(
    x = as.matrix(data_bin_NA),
    model = "lda",
    posit_y = 1,
    posit_x = c(2, 3, 4),
    w = data_bin_NA[, 6]
  )
  testthat::expect_true(cor(cbind(pred_lda[index_NA], data_bin[index_NA, 1]))[1, 2] > 0.5)
})

test_that("VIF", {
  vifs <- VIF(
    x = as.matrix(data_con_NA) * 5,
    posit_y = 1,
    posit_x = c(2, 3, 4)
  )

  collin_x <- function(y, x) {
    full_rows <- !is.na(y) & complete.cases(x)
    XXinv <- solve(crossprod(scale(x[full_rows, ], scale = F)))
    vifs <- sapply(1:ncol(XXinv), function(i) det(XXinv[i, i, drop = FALSE]) * det(XXinv[-i, -i, drop = FALSE]) / det(XXinv))
    vifs2 <- sapply(1:ncol(XXinv), function(i) Reduce("*", svd(XXinv[i, i, drop = FALSE])$d) * Reduce("*", svd(XXinv[-i, -i, drop = FALSE])$d) / Reduce("*", svd(XXinv)$d))
    list(vifs, vifs2)
  }
  cc <- collin_x(data_con_NA[, 1] * 5, data_con_NA[, c(2, 3, 4)] * 5)

  testthat::expect_true(all(cor(cbind(vifs, cc[[1]], cc[[2]]))[1, ] > 0.95))
})

data <- cbind(y_true = data_con[, 1], data_con_NA, Intercept = 1, index = 1:nrow(data_con))

data_df <- data.frame(data)
data_df$y_chac <- as.character(round(pnorm(data_df$y) * 10))
data_df$y_fac <- as.factor(round(pnorm(data_df$y) * 10))
data_df$x2 <- as.factor(round(pnorm(data_df$x2) * 5))
data_df$x3 <- as.character(round(pnorm(data_df$x3) * 5))

testthat::test_that("fill_NA error", {
  expect_error(fill_NA(
    x = data_df,
    model = "lm_pred",
    posit_y = 111,
    posit_x = c("Intercept", "x2", "x3", "x4"),
    w = data_df[["weights"]]
  ))

  expect_error(fill_NA(
    x = data_df,
    model = "lm_bayes",
    posit_y = 1,
    posit_x = 111,
    w = data_df[["weights"]]
  ))

  expect_error(fill_NA(
    x = data_df,
    model = "NOTEXISTS",
    posit_y = 1,
    posit_x = c("Intercept", "x2", "x3", "x4"),
    w = data_df[["weights"]]
  ), "'arg' should be one of")

  expect_error(fill_NA(
    x = data_df,
    model = "lda",
    posit_y = 1,
    posit_x = c("Intercepty", "NOTEXITS"),
    w = data_df[["weights"]]
  ))

  expect_error(fill_NA(
    x = 1:10,
    model = "lda",
    posit_y = 1,
    posit_x = c(3, 4),
    w = data_df[["weights"]]
  ), "wrong data type")

  expect_error(fill_NA(
    x = data_df,
    model = "pmm",
    posit_y = 1,
    posit_x = c(1, 3),
    w = data_df[["weights"]]
  ), "the same variable is dependent and indepentent")

  expect_error(fill_NA(
    x = data_df,
    model = "lda",
    logreg = TRUE,
    posit_y = 1,
    posit_x = 2,
    w = data_df[["weights"]]
  ), "logreg works only")
})

testthat::test_that("fill_NA_N error", {
  expect_error(fill_NA_N(
    x = data_df,
    model = "lm_bayes",
    posit_y = 111,
    posit_x = c("Intercept", "x2", "x3", "x4"),
    w = data_df[["weights"]]
  ))

  expect_error(fill_NA_N(
    x = data_df,
    model = "lm_bayes",
    posit_y = 1,
    posit_x = 111,
    w = data_df[["weights"]]
  ))

  expect_error(fill_NA_N(
    x = data_df,
    model = "NOTEXISTS",
    posit_y = 1,
    posit_x = c("Intercept", "x2", "x3", "x4"),
    w = data_df[["weights"]]
  ), "'arg' should be one of")

  expect_error(fill_NA_N(
    x = data_df,
    model = "lm_bayes",
    posit_y = 1,
    posit_x = c("Intercepty", "NOTEXITS"),
    w = data_df[["weights"]]
  ), "posit_x is empty")

  expect_error(fill_NA_N(
    x = 1:10,
    model = "pmm",
    posit_y = 1,
    posit_x = c(3, 4),
    w = data_df[["weights"]]
  ), "wrong data type")

  expect_error(fill_NA_N(
    x = data_df,
    model = "pmm",
    posit_y = 1,
    posit_x = c(1, 3),
    w = data_df[["weights"]]
  ), "the same variable is dependent and indepentent")

  expect_error(fill_NA_N(
    x = data_df,
    model = "pmm",
    logreg = TRUE,
    posit_y = 12,
    posit_x = 3,
    w = data_df[["weights"]]
  ), "logreg works only")
})

test_that("dplyr pipe pipeline fill_NA_N, fill_NA accuracy - con", {
  data_df2 <- data_df %>%
    group_by(group) %>%
    do(mutate(., y_imp = fill_NA(
      x = .,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"),
      w = .[["weights"]]
    ))) %>%
    do(mutate(., y_imp2 := fill_NA_N(
      x = .,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 100
    ))) %>%
    do(mutate(., y_imp3 := fill_NA(
      x = .,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept"), w = .[["weights"]]
    ))) %>%
    do(mutate(., y_imp4 := fill_NA(
      x = .,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept")
    ))) %>%
    do(mutate(., y_imp5 := fill_NA_N(
      x = .,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 10
    ))) %>%
    do(mutate(., y_imp6 := fill_NA_N(
      x = .,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 100
    ))) %>%
    ungroup() %>%
    mutate(y_imp7 := fill_NA_N(
      x = .,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    )) %>%
    group_by(group) %>%
    do(mutate(., y_imp8 = fill_NA_N(
      x = .,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 100
    ))) %>%
    ungroup() %>%
    mutate(y_imp9 = fill_NA(
      x = .,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    )) %>%
    mutate(y_imp10 = fill_NA(
      x = .,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    )) %>%
    mutate(y_imp11 = fill_NA(
      x = .,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4")
    )) %>%
    mutate(y_imp12 = fill_NA(
      x = .,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4")
    )) %>%
    mutate(y_imp13 = fill_NA_N(
      x = .,
      model = "pmm",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    )) %>%
    group_by(group) %>%
    do(mutate(., y_imp14 = fill_NA_N(
      x = .,
      model = "pmm",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    ))) %>%
    ungroup() %>%
    mutate(y_imp15 = fill_NA_N(
      x = .,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 100
    )) %>%
    mutate(y_imp16 = fill_NA_N(
      x = .,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 100
    ))

  testthat::expect_true(all(cor(data_df2[index_NA, c("y_true", colnames(data_df2)[grep("y_imp", colnames(data_df2))])])[, 1] > 0.3))
})


data_DT <- as.data.table(data_df)

test_that("data.table pipeline fill_NA_N, fill_NA accuracy - con", {
  data_DT[, y_imp := fill_NA(
    x = .SD,
    model = "lm_pred",
    posit_y = "y",
    posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
  ), by = .(group)] %>%
    .[, y_imp2 := fill_NA_N(
      x = .SD,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]], k = 100
    ), by = .(group)] %>%
    .[, y_imp3 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept"), w = .SD[["weights"]]
    ), by = .(group)] %>%
    .[, y_imp4 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept")
    ), by = .(group)] %>%
    .[, y_imp5 := fill_NA_N(
      x = .SD,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]], k = 10
    ), by = .(group)] %>%
    .[, y_imp6 := fill_NA_N(
      x = .SD,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 100
    ), by = .(group)] %>%
    .[, y_imp7 := fill_NA_N(
      x = .SD,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    )] %>%
    .[, y_imp8 := fill_NA_N(
      x = .SD,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]], k = 100
    ), by = .(group)] %>%
    .[, y_imp9 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )] %>%
    .[, y_imp10 := fill_NA(
      x = .SD,
      model = "lm_bayes",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )] %>%
    .[, y_imp11 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4")
    )] %>%
    .[, y_imp12 := fill_NA(
      x = .SD,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4")
    )] %>%
    .[, y_imp13 := fill_NA_N(
      x = .SD,
      model = "pmm",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    )] %>%
    .[, y_imp14 := fill_NA_N(
      x = .SD,
      model = "pmm",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 10
    ), by = .(group)] %>%
    .[, y_imp15 := fill_NA_N(
      x = .SD,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), k = 100
    )] %>%
    .[, y_imp16 := fill_NA_N(
      x = .SD,
      model = "lm_noise",
      posit_y = "y",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]], k = 100
    )]

  testthat::expect_true(all(cor(data_DT[index_NA, c("y_true", colnames(data_DT)[grep("y_imp", colnames(data_DT))]), with = FALSE])[, 1] > 0.3))
})


data <- cbind(y_true = data_disc[, 1], data_disc_NA, Intercept = 1, index = 1:nrow(data_disc))

data_df <- data.frame(data)
data_df$y_chac <- as.character(data_df$y)
data_df$y_fac <- as.factor(data_df$y)
data_df$y_logi <- as.logical(data_df$y >= 5)
data_df$y_true_bin <- as.logical(data_df$y_true >= 5)
data_df$x2 <- as.factor(round(pnorm(data_df$x2) * 3))
data_df$x3 <- as.character(round(pnorm(data_df$x3) * 3))
data_df$x4_bhalf <- data_df$x4 > 0.5

testthat::test_that("dplyr pipe pipeline fill_NA_N, fill_NA - dis", {
  data_df2 <- data_df %>%
    mutate(y_imp1 = fill_NA_N(
      x = .,
      model = "lm_noise",
      posit_y = "y_fac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 100
    )) %>%
    mutate(y_imp2 = fill_NA_N(
      x = .,
      model = "lm_noise",
      posit_y = "y_chac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]], k = 100
    )) %>%
    mutate(y_imp3 = fill_NA(
      x = .,
      model = "lda",
      posit_y = "y_chac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    )) %>%
    mutate(y_imp4 = fill_NA(
      x = .,
      model = "lda",
      posit_y = "y_fac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    ))  %>%
    mutate(y_imp5 = fill_NA(
      x = .,
      model = "lm_noise",
      posit_y = "y_chac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    )) %>%
    mutate(y_imp6 = fill_NA(
      x = .,
      model = "lm_noise",
      posit_y = "y_fac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .[["weights"]]
    ))

  dd <- data_df2[index_NA, colnames(data_df2)[grep("y_imp", colnames(data_df2))]]
  dd_y <- data_df2[index_NA, "y_true"]

  # Better than naive
  testthat::expect_true(all(lapply(dd, function(x) 100 * mean(x == dd_y))  > 10))

})

data_DT <- as.data.table(data_df)

testthat::test_that("data.table pipeline fill_NA/fill_NA_N accuracy - disc", {
  data_DT[, y_imp := fill_NA(
    x = .SD,
    model = "lda",
    posit_y = "y_fac",
    posit_x = c("x2", "x3", "x4", "x5")
  ), by = .(group)] %>%
    .[, y_imp2 := fill_NA(
      x = .SD,
      model = "lda",
      posit_y = "y_fac",
      posit_x = c("x2", "x3", "x4", "x5")
    )] %>%
    .[, y_imp3 := fill_NA(
      x = .SD,
      model = "lda",
      posit_y = "y_chac",
      posit_x = c("x2", "x3", "x4", "x5")
    ), by = .(group)] %>%
    .[, y_imp4 := fill_NA(
      x = .SD,
      model = "lda",
      posit_y = "y_chac",
      posit_x = c("x2", "x3", "x4", "x5")
    )] %>%
    .[, y_imp5 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y_chac",
      posit_x = c("x2", "x3", "x4", "x5")
    ), by = .(group)] %>%
    .[, y_imp6 := fill_NA(
      x = .SD,
      model = "lm_pred",
      posit_y = "y_chac",
      posit_x = c("x2", "x3")
    ), by = .(group)] %>%
    .[, y_imp8 := fill_NA(
      x = .SD,
      model = "lda",
      posit_y = "y_chac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )] %>%
    .[, y_imp9 := fill_NA(
      x = .SD,
      model = "lda",
      posit_y = "y_fac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )] %>%
    .[, y_imp10 := fill_NA_N(
      x = .SD,
      model = "pmm",
      posit_y = "y_chac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )] %>%
    .[, y_imp11 := fill_NA_N(
      x = .SD,
      model = "pmm",
      posit_y = "y_fac",
      posit_x = c("Intercept", "x2", "x3", "x4"), w = .SD[["weights"]]
    )]

  # Better than naive
  testthat::expect_true(all(data_DT[index_NA, c("y_true", colnames(data_DT)[grep("y_imp", colnames(data_DT))]), with = FALSE] %>%
                              .[, lapply(.SD, function(x) 100 * mean(y_true == x))] > 10))

})
