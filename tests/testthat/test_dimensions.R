source("main.R")

context("Checking dimension of output objects.")

test_that("Dimension of snp datasets is compatible with those supplied", {
  expect_equal(dim(sim_snps$snps), c(n, p))
  expect_equal(dim(repl_snps$snps), c(n, p))
})


test_that("Dimension of phenotype datasets is compatible with those supplied", {
  expect_equal(dim(sim_phenos$phenos), c(n, d))
  expect_equal(dim(repl_phenos$phenos), c(n, d))
})
