# File: testthat.R
# Purpose: Runs the testthat tests during R CMD check.
library(testthat)
library(cfisher)

test_check("cfisher")
