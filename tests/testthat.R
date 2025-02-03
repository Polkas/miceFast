Sys.setenv("OMP_THREAD_LIMIT" = 2)

library(testthat)
library(miceFast)
library(data.table)
library(dplyr)
library(magrittr)

test_check("miceFast")
