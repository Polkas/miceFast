set.seed(1234)

data(airquality)
data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
weights <- rgamma(nrow(data), 3, 3) # a numeric vector - positive values
groups <- as.numeric(airquality[, 5]) # a numeric vector not integers - positive values
groups2 <- as.numeric(sample(1:3, nrow(data), replace = T)) # a numeric vector not integers - positive values

testthat::test_that("constructor and destructor", {
  testthat::expect_error(model <- new(miceFast), NA)
  rm(model)
  gc()
})

testthat::test_that("set_data/get_data", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  data2 <- data
  attributes(data2) <- NULL
  attr(data2, "dim") <- dim(data)
  testthat::expect_identical(model$get_data(), data2)
  rm(model)
})

testthat::test_that("set_g/get_g", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_g(groups) # providing by a reference
  testthat::expect_identical(as.vector(model$get_g()), groups)
  rm(model)
})


testthat::test_that("sort_byg", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_g(groups2) # providing by a reference
  testthat::expect_warning(model$sort_byg())
  testthat::expect_false(identical(as.vector(model$get_index()), seq_len(nrow(data))))
  testthat::expect_true(model$is_sorted_byg())
  rm(model)
})

testthat::test_that("set_w/get_w", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_w(weights) # providing by a reference
  testthat::expect_identical(as.vector(model$get_w()), weights)
  rm(model)
})

testthat::test_that("get_index", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  testthat::expect_identical(as.integer(model$get_index()), seq_len(nrow(data)))
  rm(model)
})

testthat::test_that("get_models/get_model", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  testthat::expect_identical(model$get_models(1), "lm_pred or lm_bayes or lm_noise or pmm")
  testthat::expect_identical(model$get_model(1), "lm_pred")
  testthat::expect_identical(model$get_models(3), "no NA values for the dependent variable")
  testthat::expect_identical(model$get_model(3), "no NA values for the dependent variable")
  data2 <- data
  groups2a <- groups
  data2 <- cbind(data2, groups2a)
  model$set_data(data2) # providing by a reference
  testthat::expect_identical(model$get_models(8), "no NA values for the dependent variable")
  testthat::expect_identical(model$get_model(8), "no NA values for the dependent variable")
  groups2a[1] <- NA
  data3 <- cbind(data, groups2a)
  model$set_data(data3) # providing by a reference
  testthat::expect_identical(model$get_models(8), "recommended lda or (lm_pred,lm_bayes,lm_noise, pmm - remember to round results if needed)")
  testthat::expect_identical(model$get_model(8), "lda")
  rm(model)
})

testthat::test_that("impute/impute_N - acc", {
  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  # impute adapt to provided parmaters like w or g
  # Warning - if data is not sorted increasingly by the g then it would be done automatically
  # during a first imputation
  # Simple mean - permanent imputation at the object and data

  imp_miceFast <- model$impute_N("lm_bayes", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_noise <- model$impute_N("lm_noise", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_pred <- model$impute("lm_pred", 1, c(3, 4, 5, 6))$imputations
  imp_miceFast_pmm <- model$impute_N("pmm", 1, c(3, 4, 5, 6), 1)$imputations

  testthat::expect_true(mean((imp_miceFast - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_noise - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_pmm - imp_miceFast_pred)**2) < 1000)

  rm(model)
})

testthat::test_that("impute/impute_N  - w + g - acc", {
  data(airquality)
  data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
  weights <- rgamma(nrow(data), 3, 3) # a numeric vector - positive values
  groups <- as.numeric(airquality[, 5]) # a numeric vector not integers - positive values
  groups2 <- as.numeric(sample(1:3, nrow(data), replace = T)) # a numeric vector not integers - positive values

  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_w(weights) # providing by a reference
  model$set_g(groups) # providing by a reference
  # impute adapt to provided parmaters like w or g
  # Warning - if data is not sorted increasingly by the g then it would be done automatically
  # during a first imputation
  # Simple mean - permanent imputation at the object and data

  imp_miceFast <- model$impute_N("lm_bayes", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_noise <- model$impute_N("lm_noise", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_pred <- model$impute("lm_pred", 1, c(3, 4, 5, 6))$imputations
  imp_miceFast_pmm <- model$impute_N("pmm", 1, c(3, 4, 5, 6), 1)$imputations

  testthat::expect_true(mean((imp_miceFast - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_noise - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_pmm - imp_miceFast_pred)**2) < 1000)

  rm(model)
})

testthat::test_that("impute/impute_N  - g - acc", {
  data(airquality)
  data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
  groups <- as.numeric(airquality[, 5]) # a numeric vector not integers - positive values
  groups2 <- as.numeric(sample(1:3, nrow(data), replace = T)) # a numeric vector not integers - positive values

  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_g(groups) # providing by a reference
  # impute adapt to provided parmaters like w or g
  # Warning - if data is not sorted increasingly by the g then it would be done automatically
  # during a first imputation
  # Simple mean - permanent imputation at the object and data

  imp_miceFast <- model$impute_N("lm_bayes", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_noise <- model$impute_N("lm_noise", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_pred <- model$impute("lm_pred", 1, c(3, 4, 5, 6))$imputations
  imp_miceFast_pmm <- model$impute_N("pmm", 1, c(3, 4, 5, 6), 1)$imputations

  testthat::expect_true(mean((imp_miceFast - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_noise - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_pmm - imp_miceFast_pred)**2) < 1000)

  rm(model)
})

testthat::test_that("impute/impute_N  - w - acc", {
  data(airquality)
  data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
  weights <- rgamma(nrow(data), 3, 3) # a numeric vector - positive values

  model <- new(miceFast)
  model$set_data(data) # providing by a reference
  model$set_w(weights) # providing by a reference
  # impute adapt to provided parmaters like w or g
  # Warning - if data is not sorted increasingly by the g then it would be done automatically
  # during a first imputation
  # Simple mean - permanent imputation at the object and data

  imp_miceFast <- model$impute_N("lm_bayes", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_noise <- model$impute_N("lm_noise", 1, c(3, 4, 5, 6), 10000)$imputations
  imp_miceFast_pred <- model$impute("lm_pred", 1, c(3, 4, 5, 6))$imputations
  imp_miceFast_pmm <- model$impute_N("pmm", 1, c(3, 4, 5, 6), 1)$imputations

  testthat::expect_true(mean((imp_miceFast - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_noise - imp_miceFast_pred)**2) < 1)
  testthat::expect_true(mean((imp_miceFast_pmm - imp_miceFast_pred)**2) < 1000)

  rm(model)
})

testthat::test_that("Update , impute/impute_N  - w + g", {
  data(airquality)
  data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
  weights <- rgamma(nrow(data), 3, 3) # a numeric vector - positive values
  groups2 <- as.numeric(sample(1:3, nrow(data), replace = T)) # a numeric vector not integers - positive values

  model <- new(miceFast)
  model$set_data(data) # providing data by a reference
  model$set_w(weights) # providing by a reference
  model$set_g(groups2) # providing by a reference

  # impute adapt to provided parameters like w or g
  # Simple mean - permanent imputation at the object and data
  suppressWarnings(model$update_var(1, model$impute("lm_pred", 1, c(2))$imputations))
  model$update_var(2, model$impute_N("lm_bayes", 2, c(1, 3, 4, 5, 6), 10)$imputations)

  testthat::expect_identical(model$which_updated(), c(1L, 2L))

  rm(model)
})

testthat::test_that("data order", {
  data(airquality)
  data <- cbind(as.matrix(airquality[, -5]), intercept = 1, index = 1:nrow(airquality))
  weights <- rgamma(nrow(data), 3, 3) # a numeric vector - positive values
  groups2 <- as.numeric(sample(1:3, nrow(data), replace = T)) # a numeric vector not integers - positive values
  data_origin <- cbind(data, groups2, weights)

  model <- new(miceFast)
  model$set_data(data) # providing data by a reference
  model$set_w(weights) # providing by a reference
  model$set_g(groups2) # providing by a reference

  suppressWarnings(model$sort_byg())

  data_new <- cbind(model$get_data(), model$get_g(), model$get_w())[order(model$get_index()), ]

  dims <- dim(data_origin)
  attributes(data_origin) <- NULL
  attr(data_origin, "dim") <- dims

  testthat::expect_identical(data_origin, data_new)

  rm(model)
})
