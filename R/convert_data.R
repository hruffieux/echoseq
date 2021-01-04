# This file is part of the `echoseq` R package:
#     https://github.com/hruffieux/echoseq
#

#' Convert SNP data to a list_snps object.
#'
#' This function creates a list_snps object from a SNP data matrix supplied as
#' input.
#'
#' @param snps SNP data matrix without missing values. The entries must be 0, 1
#'   or 2.
#'
#' @return An object of class "\code{list_snps}".
#'  \item{snps}{Matrix containing the original SNP data.}
#'  \item{vec_maf}{Vector containing the SNP minor allele frequencies.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 7500
#' snps <- matrix(sample(0:2, size = n * p, replace = TRUE), nrow = n)
#'
#' list_snps <- convert_snps(snps)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{convert_phenos}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
convert_snps <- function(snps) {

  check_structure_(snps, "matrix", "numeric")
  stopifnot(all(snps %in% 0:2))


  if (is.null(rownames(snps))) {
    rownames(snps) <- paste0("ind_", 1:nrow(snps))

  }

  if (is.null(colnames(snps))) {
    colnames(snps) <- paste0("snp_", 1:ncol(snps))
  }

  vec_maf <- apply(snps, 2, mean)/2

  list_snps <- create_named_list_(snps, vec_maf)
  class(list_snps) <- "list_snps"

  list_snps

}



#' Convert phenotypic data to a list_phenos object.
#'
#' This function creates a list_phenos object from a data matrix supplied as
#' input.
#'
#' @param phenos Matrix of phenotypes (rows observations, columns phenotypic
#'   variables), without missing values.
#'
#' @return An object of class "\code{list_phenos}".
#'  \item{phenos}{Matrix containing the original phenotypic data.}
#'  \item{var_err}{Vector containing the sample phenotypic variances.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; d <- 1000
#' phenos <- matrix(rnorm(n * d), nrow = n)
#'
#' list_phenos <- convert_phenos(phenos)
#'
#' @seealso \code{\link{convert_snps}}, \code{\link{generate_snps}},
#'   \code{\link{replicate_real_snps}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
convert_phenos <- function(phenos) {

  check_structure_(phenos, "matrix", "numeric")

  if (is.null(rownames(phenos))) {
    rownames(phenos) <- paste0("ind_", 1:nrow(phenos))

  }

  if (is.null(colnames(phenos))) {
    colnames(phenos) <- paste0("pheno_", 1:ncol(phenos))
  }

  var_err <- apply(phenos, 2, var)
  ind_bl <- NULL

  list_phenos <- create_named_list_(phenos, var_err, ind_bl)

  class(list_phenos) <- "list_phenos"
  list_phenos
}
