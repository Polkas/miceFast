testthat::test_that("neibo k=1", {
  testthat::expect_identical(as.vector(neibo(1:1000, 1:10, 1)), as.numeric(1:10))
  testthat::expect_identical(as.vector(neibo(c(11, 2, 33, 111, 1, 76), c(100, 10, 4, 29), 1)), c(111, 11, 2, 33))
})

testthat::test_that("neibo k=2", {
  set.seed(1234)
  testthat::expect_identical(as.vector(neibo(1:1000, 1:10, 2)), c(1, 3, 4, 5, 6, 7, 7, 8, 10, 11))
})

testthat::test_that("neibo repro", {
  set.seed(1234)
  repro1 <- neibo(1:1000, 1:10, 2)
  set.seed(1234)
  repro2 <- neibo(1:1000, 1:10, 2)

  testthat::expect_identical(repro1, repro2)
})

testthat::test_that("neibo max min", {
  set.seed(1234)
  repro1 <- neibo(1:1000, 1:10, 100000)
  testthat::expect_true(max(repro1) <= 1000)

  set.seed(1234)
  repro2 <- neibo(1:1000, 1:10, -1)
  testthat::expect_identical(as.numeric(1:10), as.vector(repro2))
})
