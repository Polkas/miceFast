collin_x <- function(y, x) {
  full_rows <- !is.na(y) & complete.cases(x)
  XXinv <- solve(crossprod(scale(x[full_rows, ], scale = F)))

  vifs <- sapply(1:ncol(XXinv), function(i) det(XXinv[i, i, drop = FALSE]) * det(XXinv[-i, -i, drop = FALSE]) / det(XXinv))

  vifs2 <- sapply(1:ncol(XXinv),
                  function(i) Reduce("*", svd(XXinv[i, i, drop = FALSE])$d) *
                    Reduce("*", svd(XXinv[-i, -i, drop = FALSE])$d) / Reduce("*", svd(XXinv)$d)
  )
  list(vifs, vifs2)
}

airquality2 <- airquality
airquality2$Temp2 <- airquality2$Temp**2
data(air_miss)

airquality_mat <- as.matrix(airquality)

test_that("VIF correct", {
  vif_R <- collin_x(airquality2[, 1], as.matrix(airquality2[, -1]))
  model <- new(miceFast)
  model$set_data(as.matrix(airquality2))
  vif_det_miceFast <- model$vifs(1, c(2, 3, 4, 5, 6, 7))
  all_vifs_sd <- apply(do.call(rbind, list(vif_R[[1]], vif_R[[2]], as.vector(vif_det_miceFast))), 2, sd)
  expect_true(all(all_vifs_sd < 0.01))
})

testthat::test_that("VIF error", {
  expect_error(VIF(airquality2, 1, NULL, correct = 2), "is.logical\\(correct\\) is not TRUE")
  expect_error(VIF(airquality2, 1, NULL, correct = FALSE), "at least two independent variables should be provided")
  expect_error(VIF(airquality2, 1, 1, correct = FALSE), "the same variable is dependent and indepentent")
  expect_error(VIF(airquality2, 2, 22:23, correct = FALSE), "posit_x %in% ")
  expect_error(VIF(airquality2, 22, 2:3, correct = FALSE), "posit_y %in%")
  expect_error(VIF(2, 22, 2:3, correct = FALSE), "wrong data type - it should be data.frame, matrix or data.table")
  # expect_error(VIF(cbind(1, airquality2), 2, c(1, 3, 4), correct = TRUE), "Do not include an intercept")
})

test_that("VIF matrix", {
  expect_true(all(VIF(airquality_mat, 1, c(2:5), correct = FALSE) >=
                    VIF(airquality_mat, "Ozone", c("Solar.R", "Wind", "Temp", "Month"),  correct = TRUE))
              )
})

test_that("VIF data.frame", {
  expect_true(all(VIF(airquality2, 1, c(2:5), correct = FALSE) >=
                    VIF(airquality2, "Ozone", c("Solar.R", "Wind", "Temp", "Month"), correct = TRUE))
             )
})

test_that("VIF data.table", {
  setDT(air_miss)
  expect_true(any(VIF(air_miss, 1, c(2:5)) >= 1))
  expect_true(any(VIF(air_miss, "Ozone", c("Solar.R", "Wind", "Temp", "Day")) >= 1))
  expect_true(all(
    VIF(air_miss, 1, c(2:5), correct = FALSE) >=
      VIF(air_miss, "Ozone", c("Solar.R", "Wind", "Temp", "Day"), correct = TRUE))
  )
})
