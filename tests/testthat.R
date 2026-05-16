# This file is part of the standard testthat infrastructure.
# It is run during R CMD check to execute all tests in tests/testthat/

# CRAN policy: a package must use at most two cores during checks.
# Cap OpenMP threads before loading choicer so its compiled code respects the limit.
# R_DATATABLE_NUM_THREADS covers data.table's separate thread pool.
Sys.setenv(OMP_THREAD_LIMIT = "2")
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(R_DATATABLE_NUM_THREADS = "2")

library(testthat)
library(choicer)

test_check("choicer")
