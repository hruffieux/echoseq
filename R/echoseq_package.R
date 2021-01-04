#' echoseq: a package for faithful replications and simulations of genetic
#' variants, molecular expression levels and other phenotypic data
#'
#' The echoseq package provides functions to emulate molecular quantitative trait
#' locus data and clinical data associated with genetic variants under
#' user-specified association patterns. The data can be the result of pure
#' simulation or can replicate real datasets supplied by the user (which may
#' replace real data when these cannots be shared for diverse privacy reasons).
#' The data generation schemes are based on generally accepted principles of
#' population genetics (Hardy--Weinberg equilibrium, linkage-disequilibrium,
#' natural selection, pleiotropic control, sparsity assumptions).
#'
#' @section echoseq functions: convert_snps, convert_phenos, generate_dependence,
#' generate_phenos, generate_snps, replicate_real_phenos, replicate_real_snps.
#'
#' @docType package
#' @name echoseq-package
#' @importFrom stats cor qnorm quantile rbeta rbinom rlnorm rnorm rpois runif
#' sd setNames var
#' @importFrom extraDistr rtpois
#' @importFrom rlang .data
#' @importFrom dplyr %>% mutate group_by
NULL
