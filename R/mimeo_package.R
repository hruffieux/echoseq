#' mimeo: a package for faithful replications and simulations of genetic 
#' variants, molecular expression levels and other phenotypic data
#'
#' The mimeo package provides functions to emulate molecular quantitative trait
#' locus data and clinical data associated with genetic variants under
#' user-specified association patterns. The data can be the result of pure
#' simulation or can replicate real datasets supplied by the user (which may
#' replace real data when these cannot be shared for diverse privacy reasons).
#' The data generation schemes are based on generally accepted principles of
#' population genetics (Hardy--Weinberg equilibrium, linkage-disequilibrium,
#' natural selection, pleiotropic control, sparsity assumptions).
#'
#' @section mimeo functions: generate_dependence, generate_phenos,
#' generate_snps, replicate_real_phenos, replicate_real_snps.
#'
#' @docType package
#' @name mimeo-package
#' @importFrom stats cor qnorm rbeta rbinom rnorm runif setNames var
NULL
